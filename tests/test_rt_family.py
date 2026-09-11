"""Tests for the retention-time family panel.

Synthetic chromatograms only. The cases here are the ones where a wrong answer
would be silent: an inversion missed, a window centred on the wrong peak, or a
family with nothing to compare reported as consistent.
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from spectrum_annotator_ddzby import (
    plot_rt_family_panel,
    RTFamilyTrace,
    choose_rt_window,
    find_inversions,
    estimate_sialic_step,
)

RT = np.linspace(40.0, 110.0, 700)


def _peak(centre, height=1.0, width=0.5):
    return height * np.exp(-0.5 * ((RT - centre) / width) ** 2)


def _ordered_family():
    """Correctly ordered: each extra sialic acid elutes ~9 min later."""
    return [
        RTFamilyTrace('N4H5', RT, _peak(60.0), n_sialic=0),
        RTFamilyTrace('N4H5A1', RT, _peak(69.0), n_sialic=1, is_psm=True),
        RTFamilyTrace('N4H5A2', RT, _peak(78.0), n_sialic=2),
    ]


class TestInversions:
    def test_correctly_ordered_family_has_none(self):
        assert find_inversions(_ordered_family()) == []

    def test_a_sialylated_form_eluting_early_is_caught(self):
        bad = [
            RTFamilyTrace('N4H5', RT, _peak(70.0), n_sialic=0),
            RTFamilyTrace('N4H5A1', RT, _peak(60.0), n_sialic=1),
        ]
        assert find_inversions(bad) == [('N4H5A1', 'N4H5')]

    def test_is_an_ordering_test_with_no_tolerance(self):
        """One second of inversion is still an inversion. The criterion is
        ordering, so inventing a tolerance would weaken it into the outlier
        test it is meant to be independent of."""
        bad = [
            RTFamilyTrace('N4H5', RT, _peak(70.0), n_sialic=0),
            RTFamilyTrace('N4H5A1', RT, _peak(69.9), n_sialic=1),
        ]
        assert len(find_inversions(bad)) == 1

    def test_traces_without_a_peak_are_skipped_not_assumed_consistent(self):
        traces = _ordered_family() + [
            RTFamilyTrace('N4H5A3', RT, np.zeros_like(RT), n_sialic=3)]
        assert find_inversions(traces) == []
        assert traces[-1].apex_rt is None


class TestSialicStep:
    def test_measures_the_step(self):
        step = estimate_sialic_step(_ordered_family())
        assert step == pytest.approx(9.0, abs=0.3)

    def test_single_level_family_returns_none_not_zero(self):
        """Nothing to compare. None says so; 0.0 would claim a measured step."""
        one = [RTFamilyTrace('N4H5', RT, _peak(60.0), n_sialic=0),
               RTFamilyTrace('N4H5F1', RT, _peak(61.0), n_sialic=0)]
        assert estimate_sialic_step(one) is None


class TestWindow:
    def test_centres_on_the_identification_not_the_biggest_peak(self):
        """Over a long gradient the same m/z often elutes more than once, and
        the global apex is frequently not the peak that produced the PSM."""
        traces = [RTFamilyTrace('N4H5A1', RT, _peak(60.0, height=1.0)
                                + _peak(95.0, height=50.0), n_sialic=1)]
        lo, hi = choose_rt_window(60.0, traces, pad_min=5.0)
        assert (lo, hi) == pytest.approx((55.0, 65.0))
        assert not (lo <= 95.0 <= hi)

    def test_falls_back_to_the_family_span_without_a_psm_time(self):
        lo, hi = choose_rt_window(None, _ordered_family(), pad_min=5.0)
        # abs=0.1 is the RT grid spacing: an apex can only land on a sampled
        # point, so 60.0 is reported as the nearest grid value.
        assert lo == pytest.approx(55.0, abs=0.1)
        assert hi == pytest.approx(83.0, abs=0.1)
        # the window spans the family, not one member
        assert lo < 60.0 and hi > 78.0

    def test_no_peaks_and_no_psm_time_gives_none(self):
        empty = [RTFamilyTrace('x', RT, np.zeros_like(RT), n_sialic=0)]
        assert choose_rt_window(None, empty) is None


class TestPanel:
    def test_standalone_creates_its_own_figure(self):
        r = plot_rt_family_panel(_ordered_family(), psm_rt=69.0)
        assert r.fig is not None
        assert r.n_traces == 3
        assert r.psm_apex_rt == pytest.approx(69.0, abs=0.2)

    def test_embeds_into_a_caller_axes(self):
        """Same contract as plot_ms1_isolation_window, so the two panels can
        share a page."""
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        r = plot_rt_family_panel(_ordered_family(), ax=ax, psm_rt=69.0)
        assert r.fig is None
        assert r.ax is ax

    def test_reports_inversions_it_drew(self):
        bad = [RTFamilyTrace('N4H5', RT, _peak(70.0), n_sialic=0),
               RTFamilyTrace('N4H5A1', RT, _peak(60.0), n_sialic=1, is_psm=True)]
        r = plot_rt_family_panel(bad, psm_rt=60.0)
        assert r.inversions == [('N4H5A1', 'N4H5')]

    def test_empty_family_does_not_raise(self):
        r = plot_rt_family_panel([], psm_rt=60.0)
        assert r.n_traces == 0
        assert r.inversions == []
