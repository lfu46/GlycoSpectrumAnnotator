"""Retention-time family panel — the elution view of a glycopeptide backbone.

A glycoform is not judged alone. Every glycan on one peptide shares that
peptide's chromatography, so the informative picture is the whole family
overlaid: each composition's extracted ion chromatogram on one axis, ordered by
sialic acid, with the annotated identification's own apex marked.

Why sialic acid orders it. Under low-pH reversed phase the peptide dominates
retention and the glycan contributes little -- with one exception. Each
additional sialic acid moves elution measurably later, and compositions
carrying the same number co-elute in a tight cluster. Measured across the 127
peptide backbones in Urminsky et al.'s curated dataset that span two or more
sialylation levels (JACS Au 2026, doi:10.1021/jacsau.6c00875, Supporting Table
S3), the shift is a median +8.71 min per sialic acid and it holds in **127 of
127 backbones**. A glycoform eluting earlier than a less-sialylated sibling on
the same peptide is therefore inconsistent with its own assignment, and that
inversion is the single strongest post-search criterion in their hands --
carrying 20% of refutations alone and 67% jointly with precursor evidence.

The same view separates a real low-sialylated glycoform from an in-source
fragment of a higher one: the fragment shares its parent's apex exactly rather
than eluting where its own composition says it should.

Two contracts this module keeps, both inherited from :mod:`ms1_window`:

* ``ax=None`` builds its own figure; passing an ``ax`` draws into the caller's
  and returns ``fig=None``, so the panel works standalone or embedded.
* **Arrays in, not a reader.** The engine never opens files. The route layer
  extracts chromatograms -- once per run, not once per identification, since
  ``mzml_utils.extract_xics`` makes a full pass per call -- and passes them here.

The rule is specific to C18 at low pH. Porous graphitised carbon, or ion
pairing with TFA, changes elution order and it does not carry.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt


# Fraction of the measured per-sialic-acid step that still counts as "same
# cluster". Urminsky et al.'s prose calls these clusters 1-2 min, but their
# delivered data puts the p90 within-family span at 3.62 min against a median
# step of 8.71 min -- a ratio of 0.42. Expressing the window as a fraction of
# the step rather than in minutes is what lets it transfer between gradients:
# a constant tuned on their ~110 min elution range would be wrong on a 180 min
# one, while the ratio is a property of the chemistry.
DEFAULT_CLUSTER_FRACTION = 0.42


@dataclass
class RTFamilyTrace:
    """One glycoform's chromatogram within a peptide's family."""

    label: str
    rt: np.ndarray
    intensity: np.ndarray
    n_sialic: int = 0
    is_psm: bool = False

    @property
    def apex_rt(self) -> Optional[float]:
        if self.intensity is None or len(self.intensity) == 0:
            return None
        if float(np.max(self.intensity)) <= 0:
            return None
        return float(self.rt[int(np.argmax(self.intensity))])


@dataclass
class RTFamilyResult:
    """Return value of :func:`plot_rt_family_panel`.

    ``inversions`` lists ``(earlier_label, later_label)`` pairs where a
    more-sialylated glycoform eluted before a less-sialylated one on the same
    peptide -- the ordering violation, not a tolerance judgement, which is why
    it needs no threshold.
    """

    fig: Optional[plt.Figure]
    ax: plt.Axes
    n_traces: int
    apex_by_label: dict
    psm_apex_rt: Optional[float]
    inversions: List[Tuple[str, str]] = field(default_factory=list)
    sialic_step_min: Optional[float] = None


def choose_rt_window(psm_rt: Optional[float], traces: Sequence[RTFamilyTrace],
                     pad_min: float = 6.0) -> Optional[Tuple[float, float]]:
    """Retention window to draw, centred on the identification's own peak.

    Centred on ``psm_rt`` and not on the strongest peak in the trace, because
    over a long gradient the same m/z commonly elutes more than once and the
    global apex is frequently not the peak that produced the identification --
    the caveat recorded in ``OGlyco_DBA``'s ``acquisition/peak_width_profile``.
    Falling back to the family's own span when no identification time is given
    is safe; falling back to a global apex would not be.
    """
    if psm_rt is not None and np.isfinite(psm_rt):
        return (psm_rt - pad_min, psm_rt + pad_min)
    apexes = [t.apex_rt for t in traces if t.apex_rt is not None]
    if not apexes:
        return None
    return (min(apexes) - pad_min, max(apexes) + pad_min)


def find_inversions(traces: Sequence[RTFamilyTrace]) -> List[Tuple[str, str]]:
    """Glycoforms eluting earlier than a less-sialylated sibling.

    Pure ordering, so no tolerance is involved and none is invented. Traces
    without a measurable apex are skipped rather than assumed consistent.
    """
    apexed = [(t, t.apex_rt) for t in traces if t.apex_rt is not None]
    out: List[Tuple[str, str]] = []
    for t, rt in apexed:
        for other, other_rt in apexed:
            if other.n_sialic < t.n_sialic and rt < other_rt:
                out.append((t.label, other.label))
                break
    return out


def estimate_sialic_step(traces: Sequence[RTFamilyTrace]) -> Optional[float]:
    """Median retention shift per sialic acid across this family, in minutes.

    Returns None when the family spans a single sialylation level, which is the
    honest answer: with one level there is nothing to compare and no step to
    measure.
    """
    by_level: dict = {}
    for t in traces:
        if t.apex_rt is not None:
            by_level.setdefault(t.n_sialic, []).append(t.apex_rt)
    if len(by_level) < 2:
        return None
    levels = sorted(by_level)
    medians = [float(np.median(by_level[k])) for k in levels]
    return float(np.polyfit(np.array(levels, dtype=float), np.array(medians), 1)[0])


def plot_rt_family_panel(
    traces: Sequence[RTFamilyTrace],
    *,
    ax=None,
    psm_rt: Optional[float] = None,
    rt_window: Optional[Tuple[float, float]] = None,
    pad_min: float = 6.0,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (7, 3),
    show_legend: bool = True,
) -> RTFamilyResult:
    """Draw one peptide backbone's glycoform chromatograms on a single axis.

    Parameters
    ----------
    traces
        One :class:`RTFamilyTrace` per glycoform. Extract them together, in a
        single pass over the run.
    ax
        Draw into this axes to embed the panel; ``None`` creates its own figure.
    psm_rt
        Retention time of the identification being reviewed. Marked, and used
        to centre the window -- see :func:`choose_rt_window`.
    rt_window
        Explicit ``(min, max)``; otherwise chosen from ``psm_rt``.

    Returns
    -------
    RTFamilyResult
    """
    traces = list(traces)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = None

    if rt_window is None:
        rt_window = choose_rt_window(psm_rt, traces, pad_min=pad_min)

    ordered = sorted(traces, key=lambda t: (t.n_sialic, t.label))
    levels = sorted({t.n_sialic for t in ordered})
    # Colour carries the sialylation level, which is the ordering variable --
    # not the trace index, which would carry nothing.
    cmap = plt.cm.viridis(np.linspace(0.08, 0.88, max(len(levels), 1)))
    level_color = {lv: cmap[i] for i, lv in enumerate(levels)}

    apex_by_label = {}
    drew = 0
    for t in ordered:
        if t.rt is None or len(t.rt) == 0:
            continue
        colour = level_color.get(t.n_sialic, '#666666')
        ax.plot(t.rt, t.intensity, lw=1.6 if t.is_psm else 1.0, color=colour,
                alpha=1.0 if t.is_psm else 0.75,
                label=f"{t.label} (S{t.n_sialic})")
        apex = t.apex_rt
        apex_by_label[t.label] = apex
        if apex is not None:
            ax.axvline(apex, ls=':', lw=0.7, color=colour, alpha=0.55)
        drew += 1

    psm_apex = None
    for t in ordered:
        if t.is_psm:
            psm_apex = t.apex_rt
            break
    if psm_rt is not None and np.isfinite(psm_rt):
        ax.axvline(psm_rt, color='#D55E00', lw=1.2, alpha=0.9)
        ax.annotate('PSM', xy=(psm_rt, ax.get_ylim()[1]), xytext=(2, -2),
                    textcoords='offset points', fontsize=7, color='#D55E00',
                    ha='left', va='top')

    inversions = find_inversions(ordered)
    step = estimate_sialic_step(ordered)
    if inversions:
        ax.text(0.02, 0.95,
                f"{len(inversions)} inversion(s): a sialylated form elutes early",
                transform=ax.transAxes, fontsize=7, color='#D55E00',
                ha='left', va='top')

    if rt_window:
        ax.set_xlim(*rt_window)
    ax.set_xlabel('Retention time (min)', fontsize=8)
    ax.set_ylabel('MS1 intensity', fontsize=8)
    ax.tick_params(labelsize=7)
    if title:
        ax.set_title(title, fontsize=8)
    if show_legend and drew:
        ax.legend(fontsize=6, frameon=False, loc='upper right')

    return RTFamilyResult(fig=fig, ax=ax, n_traces=drew,
                          apex_by_label=apex_by_label, psm_apex_rt=psm_apex,
                          inversions=inversions, sialic_step_min=step)
