"""Regression tests for the shared MS1 isolation-window primitive.

`plot_ms1_isolation_window` is the single MS1 renderer both the O-glyco review
package (`annotate_oglyco_review_package.plot_ms1_precursor`) and the N-glyco
annotator (`nglyco_annotate_spectra.render_ms1_isolation_window`) delegate to,
so its numeric behavior (ppm, MIPS offset detection, co-isolation classification)
is load-bearing for both routes. These tests use synthetic spectra only.
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from spectrum_annotator_ddzby import plot_ms1_isolation_window, MS1WindowResult


CHARGE = 3
THEO = 1234.5678
SPACING = 1.003355 / CHARGE


def _synthetic_envelope():
    """Precursor isotope envelope M+0..M+3, one co-isolation peak, background."""
    mz = [THEO + n * SPACING for n in range(4)] + [THEO - 0.5, THEO - 2.4]
    inten = [1e6, 7e5, 3e5, 1e5, 2.5e5, 5e4]
    return np.array(mz), np.array(inten)


def test_returns_result_and_figure():
    mz, inten = _synthetic_envelope()
    res = plot_ms1_isolation_window(mz, inten, THEO * (1 + 2e-6), CHARGE, 1.6,
                                    theoretical_mz=THEO)
    assert isinstance(res, MS1WindowResult)
    assert res.fig is not None          # standalone mode creates its own figure
    assert res.precursor_found is True


def test_ppm_at_monoisotopic():
    mz, inten = _synthetic_envelope()
    # precursor selected on the monoisotopic peak, +2 ppm high
    res = plot_ms1_isolation_window(mz, inten, THEO * (1 + 2e-6), CHARGE, 1.6,
                                    theoretical_mz=THEO)
    assert res.best_iso_n == 0
    assert res.best_ppm == pytest.approx(2.0, abs=0.2)


def test_mips_offset_detected():
    mz, inten = _synthetic_envelope()
    # instrument isolated the M+1 isotope instead of the monoisotopic peak
    res = plot_ms1_isolation_window(mz, inten, (THEO + SPACING) * (1 + 1e-6),
                                    CHARGE, 1.6, theoretical_mz=THEO)
    assert res.best_iso_n == 1
    assert abs(res.best_ppm) < 3.0      # small residual after MIPS correction


def test_coisolation_quantified():
    mz, inten = _synthetic_envelope()
    res = plot_ms1_isolation_window(mz, inten, THEO, CHARGE, 1.6, theoretical_mz=THEO)
    # the -0.5 Da peak sits inside a 1.6 Da window and is not a precursor isotope
    assert res.coiso_pct > 0


def test_embedded_axis_mode():
    import matplotlib.pyplot as plt
    mz, inten = _synthetic_envelope()
    fig, ax = plt.subplots()
    res = plot_ms1_isolation_window(mz, inten, THEO, CHARGE, 1.6,
                                    theoretical_mz=THEO, ax=ax)
    assert res.fig is None              # caller owns the figure when ax is supplied
    assert res.ax is ax
    plt.close(fig)


def test_theoretical_defaults_to_precursor():
    mz, inten = _synthetic_envelope()
    # no theoretical_mz -> ppm reads ~0 and MIPS is a no-op; must not raise
    res = plot_ms1_isolation_window(mz, inten, THEO, CHARGE, 1.6)
    assert res.best_iso_n == 0
    assert res.best_ppm == pytest.approx(0.0, abs=1e-6)


def test_precursor_not_detected_is_graceful():
    mz, inten = _synthetic_envelope()
    # precursor m/z far from any peak -> precursor_found False, no exception
    res = plot_ms1_isolation_window(mz, inten, THEO + 50.0, CHARGE, 1.6,
                                    theoretical_mz=THEO + 50.0)
    assert res.precursor_found is False


class TestMeasureMS1Window:
    """The headless half, added 2026-09-10.

    `measure_ms1_window` exists so per-PSM quality control can use the MS1
    numbers without building a figure. `plot_ms1_isolation_window` now
    delegates to it, so the drawn annotation and the QC column cannot drift.
    """

    def test_plotter_and_headless_agree(self):
        from spectrum_annotator_ddzby import measure_ms1_window
        mz, inten = _synthetic_envelope()
        m = measure_ms1_window(mz, inten, THEO, CHARGE, 0.7, THEO)
        r = plot_ms1_isolation_window(mz, inten, THEO, CHARGE, 0.7, THEO)
        assert r.best_iso_n == m.best_iso_n
        assert r.best_ppm == pytest.approx(m.best_ppm)
        assert r.coiso_pct == pytest.approx(m.coiso_pct if m.coiso_pct is not None else 0.0)

    def test_detects_off_by_one_monoisotope_selection(self):
        """The measurement this whole module exists for.

        An M+1 selection is not a cosmetic mass error. Several monosaccharide
        swaps sit within ~0.02 Da of a neutron, so an off-by-one selection puts
        a *different glycan composition* inside the search tolerance.
        """
        from spectrum_annotator_ddzby import measure_ms1_window
        mz, inten = _synthetic_envelope()
        m = measure_ms1_window(mz, inten, THEO + SPACING, CHARGE, 0.7, THEO)
        assert m.best_iso_n == 1
        assert abs(m.best_ppm) < 1.0

        m2 = measure_ms1_window(mz, inten, THEO + 2 * SPACING, CHARGE, 0.7, THEO)
        assert m2.best_iso_n == 2

    def test_unmeasurable_coisolation_is_none_not_zero(self):
        """A zero and an unmeasurable co-isolation mean opposite things.

        0.0 says "clean precursor, reporter ions trustworthy". None says "no
        envelope signal found, this was never checked". Collapsing the two
        would let an unchecked PSM pass as a clean one.
        """
        from spectrum_annotator_ddzby import measure_ms1_window
        m = measure_ms1_window(
            np.array([500.0]), np.array([1.0]), THEO, CHARGE, 0.7, THEO
        )
        assert m.coiso_pct is None
        assert m.precursor_found is False
        assert m.envelope_intensity == 0.0

    def test_clean_precursor_reports_zero_not_none(self):
        from spectrum_annotator_ddzby import measure_ms1_window
        mz = np.array([THEO + n * SPACING for n in range(4)])
        inten = np.array([100.0, 90.0, 50.0, 20.0])
        m = measure_ms1_window(mz, inten, THEO, CHARGE, 0.7, THEO)
        assert m.coiso_pct == pytest.approx(0.0)
        assert m.n_coiso_peaks == 0
        assert m.precursor_found is True

    def test_coisolation_ratio_is_window_restricted(self):
        """Both numerator and denominator are limited to the isolation window.

        The question being answered is "of the ions the quadrupole actually let
        through, what fraction is not my precursor" -- which is what determines
        reporter-ion contamination. Counting envelope peaks that fell outside
        the window into the denominator understates co-isolation, and
        understating is the dangerous direction.
        """
        from spectrum_annotator_ddzby import measure_ms1_window
        # M+0 and M+1 fall inside a 0.7 Da window at charge 3; M+2 and M+3 do not.
        mz = np.array([THEO + n * SPACING for n in range(4)] + [THEO + 0.25])
        inten = np.array([100.0, 90.0, 50.0, 20.0, 19.0])
        m = measure_ms1_window(mz, inten, THEO, CHARGE, 0.7, THEO)
        assert m.envelope_intensity == pytest.approx(190.0)   # not 260.0
        assert m.coiso_intensity == pytest.approx(19.0)
        assert m.coiso_pct == pytest.approx(19.0 / 190.0 * 100.0)
