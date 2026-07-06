"""MS1 isolation-window rendering — shared by the O-glyco and N-glyco annotators.

Draws the MS1 precursor region zoomed to the isolation window, classifying peaks
as precursor-envelope / co-isolation / background, with monoisotopic (MIPS)
offset detection and precursor ppm annotation.

Reader-agnostic by design: callers pass raw m/z + intensity arrays (pulled from
any ``mzml_utils.open_spectra`` reader), so this module stays inside the
annotation engine with no dependency on the reader layer. This is the single
implementation used by both:

  * O-glyco review package — ``annotate_impa_buffer_spectra.plot_ms1_precursor``
    (panel set: MS1 isolation window + HCD + EThcD)
  * N-glyco annotator — ``nglyco_annotate_spectra.render_ms1_isolation_window``
    (panel set: MS1 isolation window + HCD; N-glyco is HCD-only, no EThcD)

Keeping one renderer here means tolerance/co-isolation/MIPS logic can never drift
between the two glycopeptide types.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt


# 13C-12C mass difference; per-charge isotope spacing is this divided by charge.
ISOTOPE_SPACING_DA = 1.003355


@dataclass
class MS1WindowResult:
    """Return value of :func:`plot_ms1_isolation_window`.

    ``fig`` is None when the caller supplied its own ``ax`` (embedded panel);
    the caller then owns the enclosing figure.
    """

    fig: Optional[plt.Figure]
    ax: plt.Axes
    best_ppm: float          # precursor ppm at the best-matching isotope offset
    best_iso_n: int          # MIPS offset (0 = monoisotopic peak was selected)
    coiso_pct: float         # co-isolation intensity as % of the precursor envelope
    precursor_found: bool


def plot_ms1_isolation_window(
    ms1_mz,
    ms1_intensity,
    precursor_mz,
    charge,
    isolation_width,
    theoretical_mz=None,
    *,
    ax=None,
    ms1_scan=None,
    trigger_label='HCD',
    figsize=(7, 3),
):
    """Render the MS1 isolation window around a precursor.

    Parameters
    ----------
    ms1_mz, ms1_intensity : array-like
        Full MS1 spectrum peaks (e.g. ``reader.get_spectrum(ms1_scan).mz`` /
        ``.intensity``).
    precursor_mz : float
        Isolation-center m/z the instrument selected.
    charge : int
        Precursor charge.
    isolation_width : float
        Isolation window width in Da (drawn as the shaded band).
    theoretical_mz : float, optional
        Theoretical monoisotopic precursor m/z. For MSFragger data use
        ``(Calculated Peptide Mass + z * proton) / z`` (the full glycopeptide
        mass, MIPS-corrected); for OPair data use ``FragmentCalculator.precursor_mz``.
        If None, defaults to ``precursor_mz`` so ppm reads ~0 and MIPS detection
        is a no-op (the window still renders).
    ax : matplotlib Axes, optional
        Draw onto this axes to embed the panel in a multi-panel figure. If None,
        a standalone ``(fig, ax)`` is created and returned in the result.
    ms1_scan : int, optional
        MS1 scan number, used only in the title.
    trigger_label : str
        Right-hand side of the ``"MS1 scan N -> {trigger_label}"`` title, e.g.
        ``'HCD scan 13399'`` or ``'HCD 13399 + EThcD 13401-13403'``.

    Returns
    -------
    MS1WindowResult
    """
    mz = np.asarray(ms1_mz, dtype=float)
    intensity = np.asarray(ms1_intensity, dtype=float)
    if theoretical_mz is None:
        theoretical_mz = precursor_mz
    isotope_spacing = ISOTOPE_SPACING_DA / charge

    # Determine best isotope match between selected and theoretical (MIPS check:
    # detects when M+1..M+3 was isolated instead of the monoisotopic peak).
    best_iso_n = 0
    best_ppm = (precursor_mz - theoretical_mz) / theoretical_mz * 1e6
    for iso_n in range(4):
        corrected_mz = theoretical_mz + iso_n * isotope_spacing
        ppm = (precursor_mz - corrected_mz) / corrected_mz * 1e6
        if abs(ppm) < abs(best_ppm):
            best_ppm = ppm
            best_iso_n = iso_n

    half_iso = isolation_width / 2
    margin = 3.0
    lo = precursor_mz - half_iso - margin
    hi = precursor_mz + half_iso + margin

    mask = (mz >= lo) & (mz <= hi)
    mz_zoom = mz[mask]
    int_zoom = intensity[mask]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = None

    is_precursor = np.zeros(len(mz_zoom), dtype=bool)
    is_coisolation = np.zeros(len(mz_zoom), dtype=bool)
    coiso_pct = 0.0
    if len(mz_zoom) > 0:
        # Classify peaks: precursor envelope vs co-isolation vs background
        iso_window_lo = precursor_mz - half_iso
        iso_window_hi = precursor_mz + half_iso
        in_window = (mz_zoom >= iso_window_lo) & (mz_zoom <= iso_window_hi)
        for n in range(-1, 6):
            target = theoretical_mz + n * isotope_spacing
            is_precursor |= np.abs(mz_zoom - target) < 0.02
        is_coisolation = in_window & ~is_precursor
        is_background = ~in_window

        if is_background.any():
            ax.vlines(mz_zoom[is_background], 0, int_zoom[is_background],
                      colors='#888888', linewidth=0.6, alpha=0.6)
        if is_precursor.any():
            ax.vlines(mz_zoom[is_precursor], 0, int_zoom[is_precursor],
                      colors='#0d75bc', linewidth=1.2)
        if is_coisolation.any():
            ax.vlines(mz_zoom[is_coisolation], 0, int_zoom[is_coisolation],
                      colors='#e67e22', linewidth=1.0)

    ax.axvspan(precursor_mz - half_iso, precursor_mz + half_iso,
               alpha=0.08, color='blue')
    ax.axvline(precursor_mz, color='red', linestyle='--', linewidth=1.0, alpha=0.7)

    # Mark theoretical isotope peaks
    for iso_n in range(4):
        iso_mz = theoretical_mz + iso_n * isotope_spacing
        ax.axvline(iso_mz, color='green', linestyle=':', linewidth=0.7, alpha=0.5)

    precursor_found = False
    if len(mz_zoom) > 0:
        prec_mask = np.abs(mz_zoom - precursor_mz) < 0.02
        if prec_mask.any():
            precursor_found = True
            prec_idx = np.argmax(int_zoom * prec_mask)
            prec_int = int_zoom[prec_idx]

            iso_note = f'\nMIPS: M+{best_iso_n}' if best_iso_n > 0 else ''

            ax.annotate(
                f'Precursor\n{precursor_mz:.4f} m/z\n{best_ppm:+.1f} ppm\nz={charge}{iso_note}',
                xy=(mz_zoom[prec_idx], prec_int),
                xytext=(0.75, 0.85), textcoords='axes fraction',
                fontsize=8, color='#0d75bc', ha='center',
                arrowprops=dict(arrowstyle='->', color='#0d75bc', lw=1.0))

            # Label isotope peaks
            for iso_n in range(4):
                iso_mz = theoretical_mz + iso_n * isotope_spacing
                iso_mask = np.abs(mz_zoom - iso_mz) < 0.02
                if iso_mask.any():
                    iso_idx = np.argmax(int_zoom * iso_mask)
                    ax.text(mz_zoom[iso_idx], int_zoom[iso_idx],
                            f' M+{iso_n}', fontsize=6, va='bottom', color='#0d75bc')

            # Label co-isolation peaks
            if is_coisolation.any():
                n_coiso = int(is_coisolation.sum())
                coiso_max_idx = np.argmax(int_zoom * is_coisolation)
                coiso_pct = (int_zoom[is_coisolation].sum() / int_zoom[is_precursor].sum() * 100
                             if is_precursor.any() else 0.0)
                ax.annotate(f'Co-isolation\n{n_coiso} peaks ({coiso_pct:.0f}% RI)',
                            xy=(mz_zoom[coiso_max_idx], int_zoom[coiso_max_idx]),
                            xytext=(0.20, 0.85), textcoords='axes fraction',
                            fontsize=7, color='#e67e22', ha='center',
                            arrowprops=dict(arrowstyle='->', color='#e67e22', lw=0.8))
        else:
            ax.text(0.5, 0.5, 'Precursor peak not detected',
                    transform=ax.transAxes, ha='center', fontsize=10, color='red')

    ax.set_xlim(lo, hi)
    ax.set_xlabel('m/z', fontsize=10)
    ax.set_ylabel('Intensity', fontsize=10)
    iso_label = f' (MIPS M+{best_iso_n})' if best_iso_n > 0 else ''
    scan_prefix = f'MS1 scan {ms1_scan}' if ms1_scan is not None else 'MS1'
    ax.set_title(f'{scan_prefix} → {trigger_label}  |  '
                 f'Isolation: {isolation_width:.2f} Da  |  '
                 f'Precursor ppm: {best_ppm:+.2f}{iso_label}',
                 fontsize=9)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

    if fig is not None:
        fig.tight_layout()

    return MS1WindowResult(
        fig=fig,
        ax=ax,
        best_ppm=best_ppm,
        best_iso_n=best_iso_n,
        coiso_pct=coiso_pct,
        precursor_found=precursor_found,
    )
