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
