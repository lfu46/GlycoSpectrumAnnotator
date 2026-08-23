#!/usr/bin/env python3
"""
Spectrum Annotator for Glycopeptides

This module creates publication-quality annotated EThcD spectra for
glycopeptide identification, similar to IPSA output.

Author: Longping Fu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from typing import List, Dict, Optional, Tuple, Union
import os
import json
import re

# Import the fragment calculator (relative import for package)
from mzml_utils import (
    MatchContext as _MatchContext,
    canonical_match as _canonical_match,
    from_arrays as _match_ctx_from_arrays,
)
from .fragment_calculator import (
    FragmentCalculator,
    TheoreticalIon,
    MatchedIon,
    FalseMatchRate,
    match_peaks,
    parse_modifications_from_string,
    calculate_false_match_rate,
    calculate_annotation_statistics,
    OXONIUM_IONS
)

# =============================================================================
# COLOR SCHEME (IPSA-style colors)
# =============================================================================

# Ion type colors (IPSA palette)
ION_COLORS = {
    'b': '#0d75bc',      # Blue (b ions - HCD N-terminal)
    'y': '#be202d',      # Red (y ions - HCD C-terminal)
    'c': '#07a14a',      # Green (c ions - ETD N-terminal)
    'z': '#f79420',      # Orange (z ions - ETD C-terminal)
    'Y': '#17A589',      # Teal (Y-ladder - partial glycan retained)
    'Y0': '#9B59B6',     # Purple (Y0 - bare peptide, complete glycan loss)
    'Y1': '#17A589',     # Teal (Y-ladder - partial glycan retained)
    'oxonium': '#7E6148', # Brown (oxonium/glycan diagnostic)
    'precursor': '#666666',  # Dark gray (precursor)
    'unassigned': '#a6a6a6', # Gray (unmatched - IPSA style)
}

# Neutral loss suffix colors (slightly lighter/transparent versions)
ION_COLORS_NL = {
    'b_NL': '#0d75bc80',
    'y_NL': '#be202d80',
    'c_NL': '#07a14a80',
    'z_NL': '#f7942080',
}

# =============================================================================
# SPECTRUM ANNOTATOR CLASS
# =============================================================================

def resolve_peak_labels(labels,
                        max_labels: int = 50,
                        intensity_threshold_pct: float = 1.0,
                        mz_window: float = 15.0,
                        intensity_window_pct: float = 8.0):
    """The spectrum panel's peak-label placement, extracted from ``plot()`` 2026-08-23.

    Greedy and intensity-descending: the strongest peak is placed first and overlap-blocks
    weaker neighbours, never the other way around. A label inside another's exclusion window
    (``mz_window`` m/z x ``intensity_window_pct`` %-RI) is DROPPED, not moved; labels longer
    than 4 characters rotate 90 degrees. This is the whole overlap policy of the annotated
    spectrum -- public so a caller rendering a spectrum through its own plotting stack (a
    journal-styled panel, for instance) can let THIS algorithm choose which labels survive
    instead of re-inventing a collision rule.

    Args:
        labels: iterable of ``(mz, rel_intensity_pct, text)``.
        max_labels / intensity_threshold_pct: the same budget and floor ``plot()`` takes.
        mz_window / intensity_window_pct: the exclusion rectangle (defaults are the deployed
            15 m/z x 8 %).

    Returns:
        ``[(mz, rel_intensity_pct, text, rotation_deg), ...]`` for the survivors, in placement
        order (strongest first).
    """
    kept = []
    placed = []
    for mz, rel, text in sorted(labels, key=lambda l: -l[1]):
        if rel < intensity_threshold_pct or len(kept) >= max_labels:
            continue
        if any(abs(mz - px) < mz_window and abs(rel - py) < intensity_window_pct
               for px, py in placed):
            continue
        kept.append((float(mz), float(rel), text, 90 if len(text) > 4 else 0))
        placed.append((mz, rel))
    return kept


class SpectrumAnnotator:
    """
    Creates annotated spectrum plots for glycopeptides.
    """

    def __init__(self,
                 peptide: str,
                 modifications: List[Dict],
                 precursor_charge: int,
                 precursor_mz: float,
                 exp_mz: np.ndarray,
                 exp_intensity: np.ndarray,
                 tolerance_ppm: float = 20.0,
                 site_index: str = "",
                 gene: str = "",
                 activation_type: str = "HCD",
                 glycan_labels: Optional[Dict[int, str]] = None,
                 do_deisotope: bool = True,
                 noise_cache: Optional[Dict] = None,
                 scan_num: Optional[int] = None,
                 sn_threshold: float = 5.0,
                 confidence_level: Optional[str] = None,
                 source_file: Optional[str] = None,
                 extended_y_series: bool = True,
                 n_sequon_positions: Optional[List[int]] = None,
                 title_override: Optional[str] = None,
                 max_fragment_charge: Optional[int] = None):
        """
        Initialize the spectrum annotator.

        Args:
            peptide: Peptide sequence
            modifications: List of modification dicts
            precursor_charge: Precursor charge state
            precursor_mz: Experimental precursor m/z
            exp_mz: Experimental m/z array
            exp_intensity: Experimental intensity array
            tolerance_ppm: Mass tolerance for matching (default 20 ppm)
            site_index: Site identifier (e.g., "Q96KR1_S195")
            gene: Gene name
            activation_type: Fragmentation method - "HCD" (b/y ions) or "EThcD" (b/y/c/z ions)
            glycan_labels: Dict mapping 1-based glycan positions to composition strings
                          (e.g., {1: "HexNAc·Hex", 3: "HexNAc", 5: "HexNAc·Hex"}).
                          When provided, shows per-site composition below residues
                          and uses specific labels for fragment ions.
            do_deisotope: Whether to deisotope the spectrum before annotation
            noise_cache: Noise cache from load_noise_cache() for S/N calculation
            scan_num: Scan number for noise cache lookup
            sn_threshold: Minimum S/N for bond coverage counting (default 5.0)
            extended_y_series: Whether to generate extended Y-ion ladder for complex
                              glycans (default False = only Y0 + Y intact)
        """
        self.peptide = peptide
        self.modifications = modifications
        self.precursor_charge = precursor_charge
        self.precursor_mz = precursor_mz
        self.tolerance_ppm = tolerance_ppm
        self.site_index = site_index
        self.gene = gene
        self.activation_type = activation_type.upper()
        self.glycan_labels = glycan_labels or {}
        self.sn_threshold = sn_threshold
        self.scan_num = scan_num
        self.confidence_level = confidence_level
        self.source_file = source_file
        self.extended_y_series = extended_y_series
        self.n_sequon_positions = set(n_sequon_positions or [])  # 1-based positions of N in N-X-S/T
        self.title_override = title_override

        # Keep the RAW spectrum separately — used for visualization so full
        # isotope envelopes of matched ions remain visible in the plot, even
        # when deisotoping aggressively collapsed them (e.g. Y(intact) at 3+
        # whose averagine envelope sometimes gets mis-lumped with a
        # co-isolated neighbor and incorrectly removed during deisotoping).
        import numpy as _np
        self.raw_mz = _np.asarray(exp_mz)
        self.raw_intensity = _np.asarray(exp_intensity)

        # Deisotope spectrum before annotation (averagine-scored). Matching
        # is done against the deisotoped copy to avoid false matches at
        # isotope-peak m/z positions.
        if do_deisotope:
            from mzml_utils import deisotope
            self.exp_mz, self.exp_intensity = deisotope(
                exp_mz, exp_intensity, max_charge=precursor_charge
            )
        else:
            self.exp_mz = exp_mz
            self.exp_intensity = exp_intensity

        # Calculate S/N for each peak
        from .fragment_calculator import get_sn_array
        self.peak_sn = get_sn_array(
            self.exp_mz, self.exp_intensity, noise_cache, scan_num
        )

        # Calculate theoretical fragments. HCD defaults to 1+/2+ backbone
        # fragments to avoid excessive false-positive high-charge labels.
        # EThcD/ETD defaults to precursor charge - 1 because glycan-retaining
        # c/z ions are often multiply charged.
        if max_fragment_charge is None:
            if self.activation_type == "HCD":
                effective_max_fragment_charge = min(2, max(1, precursor_charge - 1))
            else:
                effective_max_fragment_charge = max(2, precursor_charge - 1)
        else:
            effective_max_fragment_charge = max(1, int(max_fragment_charge))

        # Calculate theoretical fragments
        self.calculator = FragmentCalculator(
            peptide, modifications, precursor_charge,
            max_fragment_charge=effective_max_fragment_charge,
        )

        # Get theoretical ions based on activation type
        if self.activation_type == "HCD":
            self.theoretical_ions = self.calculator.get_hcd_ions_flat(
                include_neutral_losses=True,
                neutral_loss_types=['H2O', 'NH3'],
                extended_y_series=self.extended_y_series,
            )
        else:
            self.theoretical_ions = self.calculator.get_all_ions_flat(
                include_neutral_losses=True,
                neutral_loss_types=['H2O', 'NH3'],
                extended_y_series=self.extended_y_series,
            )

        # Filter precursor ions based on activation type
        if self.activation_type == "HCD":
            self.theoretical_ions = [
                ion for ion in self.theoretical_ions
                if ion.ion_type != 'precursor'
            ]
        elif self.activation_type.upper() in ("ETHCD", "ETD"):
            self.theoretical_ions = [
                ion for ion in self.theoretical_ions
                if ion.ion_type != 'precursor'
                or ion.charge < self.precursor_charge
                or '-' in ion.annotation
            ]

        # Canonical match layer (Stage 1b Commit B) — primary + raw fallback.
        # Behaviour is identical to the previous two-call pattern:
        #   1. match_peaks(theo, exp_mz, exp_int, match_isotopes=False)
        #      (deisotoped peaks → monoisotopic only; no isotope walk).
        #   2. for ions still unmatched by (type, number, charge), a second
        #      match_peaks against the RAW spectrum with match_isotopes=False.
        # This recovers ions whose true monoisotopic peak was eaten by
        # over-aggressive deisotoping (common for large glycopeptides at 3+
        # where averagine envelopes are wide and the deisotoper mis-identifies
        # the cluster's lowest peak as the mono).
        _match_ctx = _match_ctx_from_arrays(
            exp_mz=self.exp_mz, exp_intensity=self.exp_intensity,
            raw_mz=self.raw_mz, raw_intensity=self.raw_intensity,
            tolerance_ppm=tolerance_ppm,
            deisotope_applied=True,
        )
        # Pass the annotator's local match_peaks so the returned objects are
        # the extended MatchedIon subclass (with has_modification), which
        # downstream code relies on.
        self.matched_ions = _canonical_match(
            self.theoretical_ions, _match_ctx, raw_fallback=True,
            match_fn=match_peaks,
        )

        # MIPS-aware Y-type matching: when the instrument selected a +N isotope
        # as the precursor (MIPS offset = N > 0), Y(intact) and Y-ladder ions
        # in MS2 appear at m/z shifted by N × 1.003355/charge from the true
        # monoisotopic mass. Detect MIPS from the observed precursor m/z and
        # generate MIPS-shifted theoretical Y-ions for matching. This ONLY
        # applies to Y-type ions (Y(intact), Y-intermediates, Y-ladder) —
        # backbone b/y/c/z fragments come from the intact molecule and appear
        # at true-mono positions regardless of MIPS.
        ISO_SPACING = 1.003355
        PROTON = 1.007276
        observed_prec_neutral = precursor_mz * precursor_charge - precursor_charge * PROTON
        theo_prec_neutral = self.calculator.precursor_mass
        if theo_prec_neutral > 0 and precursor_mz > 0:
            delta_neutral = observed_prec_neutral - theo_prec_neutral
            mips_iso = round(delta_neutral / ISO_SPACING)
            # Only apply if MIPS > 0 and within plausible range [1..4]
            # AND the residual (non-integer part) is within 10 ppm
            residual_ppm = abs((delta_neutral - mips_iso * ISO_SPACING)
                               / theo_prec_neutral * 1e6) if theo_prec_neutral else 99
            if 1 <= mips_iso <= 4 and residual_ppm < 10:
                self.mips_offset = mips_iso
                # Build MIPS-shifted Y-ions: clone theoretical Y-ions with
                # m/z += mips × 1.003/charge, different ion_number so dedupe
                # doesn't mask them. Mark the annotation with a MIPS tag.
                from .fragment_calculator import TheoreticalIon
                mips_y_ions = []
                for t in self.theoretical_ions:
                    if t.ion_type != 'Y' or t.charge <= 0:
                        continue
                    shifted_mz = t.mz + mips_iso * ISO_SPACING / t.charge
                    mips_y_ions.append(TheoreticalIon(
                        ion_type='Y',
                        ion_number=t.ion_number,
                        charge=t.charge,
                        mz=shifted_mz,
                        sequence=t.sequence,
                        annotation=t.annotation,  # keep same annotation
                    ))
                if mips_y_ions:
                    mips_matches = match_peaks(
                        mips_y_ions, self.raw_mz, self.raw_intensity,
                        tolerance_ppm, match_isotopes=False
                    )
                    # Dedupe against the already-matched set by BOTH exp_mz
                    # AND (ion_type, ion_number, charge) identity. The m/z-only
                    # guard missed the case where a MIPS-shifted Y-ion landed
                    # on a different peak from its already-matched true-mono
                    # twin — we'd add the same logical Y-ion twice with
                    # different exp_mz, which double-counts it in scoring and
                    # confuses the label renderer.
                    matched_mz_so_far = {m.exp_mz for m in self.matched_ions}
                    matched_id_so_far = {(m.ion_type, m.ion_number, m.charge)
                                         for m in self.matched_ions}
                    for m in mips_matches:
                        if m.exp_mz in matched_mz_so_far:
                            continue
                        if (m.ion_type, m.ion_number, m.charge) in matched_id_so_far:
                            continue
                        self.matched_ions.append(m)
                        matched_mz_so_far.add(m.exp_mz)
                        matched_id_so_far.add((m.ion_type, m.ion_number, m.charge))
            else:
                self.mips_offset = 0
        else:
            self.mips_offset = 0

        # ---------------------------------------------------------------
        # Step 1: Charge-reduced precursor exclusion (EThcD/ETD only).
        # Covers z < precursor_charge only. Reviewer-approved Phase 0
        # kept as-is; Phase 0b adds Step 2 below for original-charge
        # envelope overlap that this filter cannot catch.
        # ---------------------------------------------------------------
        # Uses the calculator's canonical _get_glycan_mass() helper (which
        # shares a single source of truth with _get_glycan_positions), instead
        # of the previous name/mass heuristic — that heuristic missed any
        # complex-glycan modification whose name didn't start with 'glycan' or
        # 'hexnac' and whose mass didn't match one of three hardcoded values,
        # so multi-site or medium-mass O-glycans (365, 406, ...) silently
        # skipped the charge-reduced exclusion step.
        if self.activation_type.upper() in ("ETHCD", "ETD"):
            from .fragment_calculator import build_charge_reduced_exclusions, filter_charge_reduced
            glycan_mass = self.calculator._get_glycan_mass()
            if glycan_mass > 0:
                cr_exclusions = build_charge_reduced_exclusions(
                    self.calculator.precursor_mass, precursor_charge, glycan_mass
                )
                self.matched_ions = filter_charge_reduced(
                    self.matched_ions, cr_exclusions
                )

        # ---------------------------------------------------------------
        # Step 2 (Phase 0b): ORIGINAL-charge precursor envelope exclusion.
        # Applies to BOTH HCD and EThcD. Reviewer identified that
        # build_charge_reduced_exclusions only covers z < precursor_charge,
        # but spurious fragment labels frequently coincide with the ORIGINAL
        # precursor isotope peaks at the ORIGINAL charge (e.g. case 7's
        # c3+N2 1+ landing on 4+ precursor M+2). This step uses the OBSERVED
        # selected precursor m/z as the envelope center (not theoretical
        # precursor_mass) so MIPS offsets and small mass errors do not drift
        # the exclusion window.
        # ---------------------------------------------------------------
        from .fragment_calculator import (
            build_precursor_envelope_exclusions,
            filter_precursor_envelope_overlap,
        )
        pe_exclusions = build_precursor_envelope_exclusions(
            precursor_mz, precursor_charge,
            n_isotopes_below=2,
            n_isotopes_above=8,
        )
        self.matched_ions, self.excluded_precursor_envelope_matches = (
            filter_precursor_envelope_overlap(
                self.matched_ions, pe_exclusions, precursor_mz, precursor_charge,
            )
        )

        # ---------------------------------------------------------------
        # Step 3: Build S/N lookup for matched ions
        # ---------------------------------------------------------------
        sn_lookup = {}
        for i, mz_val in enumerate(self.exp_mz):
            sn_lookup[mz_val] = self.peak_sn[i] if i < len(self.peak_sn) else 0
        for ion in self.matched_ions:
            ion._sn = sn_lookup.get(ion.exp_mz, 0)

        # ---------------------------------------------------------------
        # Step 4 (Phase 0b reorder, Phase 0e extension): Isotope-consistency
        # flag MUST run BEFORE peak_annotations construction so the priority
        # function can see isotope_flag on each ion.
        # ---------------------------------------------------------------
        # FLAG, not FILTER. For each eligible ion we check whether the
        # M+1 isotope partner at m/z + 1.003355/charge is present in the
        # raw spectrum, and whether its ratio vs M is within rough
        # averagine expectations. Never remove a match based on this.
        #
        # Phase 0e — extended from Y-only to Y + glycan-bearing backbone
        # (b/c/y/z with has_modification=True) to catch spurious glycan-
        # bearing c/z labels like case 16 c4+N2HAF @ 1580 and case 23
        # c1+NA° @ 810. Non-glycan backbone ions are excluded because
        # their low-mass M+1 ratios are noisy, and oxonium/precursor ions
        # are excluded as before. We additionally require expected M+1
        # ratio >= 0.3 (mass >= ~375 Da) so the averagine check has enough
        # signal to be meaningful — this preserves the reviewer rule that
        # small backbone ions are not hard-filtered by isotope absence.
        _RAW_MAX_INT = float(np.max(self.raw_intensity)) if len(self.raw_intensity) else 0.0
        for ion in self.matched_ions:
            ion.isotope_flag = None
            # Only Y ions or glycan-bearing backbone (b/c/y/z) are eligible
            if ion.ion_type == 'Y':
                pass
            elif ion.ion_type in ('b', 'c', 'y', 'z') and ion.has_modification:
                pass
            else:
                continue
            if ion.charge <= 0 or _RAW_MAX_INT <= 0 or ion.exp_intensity <= 0:
                continue
            iso_mz = ion.exp_mz + 1.003355 / ion.charge
            iso_tol_da = iso_mz * tolerance_ppm / 1e6
            mask = np.abs(self.raw_mz - iso_mz) <= iso_tol_da
            # Rough averagine expected M+1/M ratio: ~0.0008 Da^-1 of mono
            # neutral mass (Senko averagine), capped at 0.85 for very
            # large glycopeptides. Charged ion mass ~= mz * z - z * proton.
            mass_da = ion.exp_mz * max(1, ion.charge) - ion.charge * 1.007276
            expected_ratio = min(0.85, 0.0008 * mass_da)
            # Reliability threshold: below ~375 Da the expected M+1 ratio
            # is <0.3, within the noise floor of real glycopeptide spectra;
            # a missing M+1 there does not distinguish real from spurious.
            if expected_ratio < 0.3:
                continue
            if mask.any():
                iso_int = float(self.raw_intensity[mask].max())
                ratio = iso_int / ion.exp_intensity if ion.exp_intensity > 0 else 0.0
                if ratio < 0.3 * expected_ratio:
                    ion.isotope_flag = 'suspicious_low_M1'
                elif ratio > 3.0 * expected_ratio:
                    ion.isotope_flag = 'suspicious_high_M1'
                else:
                    ion.isotope_flag = 'ok'
            else:
                ion.isotope_flag = 'M1_absent'

        # ---------------------------------------------------------------
        # Step 5 (Phase 0b Patch D): peak_annotations with isotope-aware
        # priority. M+1 absence is audit-only for centroided MS2 peak lists:
        # peak picking can keep a monoisotopic centroid while dropping its
        # isotope partner. Penalize only abnormal isotope ratios when an
        # isotope partner is actually detected but inconsistent.
        # ---------------------------------------------------------------
        def _ion_priority(ion):
            flag = getattr(ion, "isotope_flag", None)
            flag_penalty = 10 if flag in {
                "suspicious_low_M1", "suspicious_high_M1"
            } else 0
            if ion.ion_type in ("Y", "Y0"):
                base = 0
            elif (
                ion.ion_type in ("b", "y", "c", "z")
                and (not ion.neutral_loss or ion.neutral_loss == "glyc_full")
            ):
                base = 1
            elif ion.ion_type in ("b", "y", "c", "z") and ion.neutral_loss:
                base = 2
            elif ion.ion_type == "oxonium":
                base = 3
            else:
                base = 4
            return base + flag_penalty

        self.matched_mz_set = set()
        self.peak_annotations = {}
        for ion in self.matched_ions:
            self.matched_mz_set.add(ion.exp_mz)
            existing = self.peak_annotations.get(ion.exp_mz)
            if existing is None or _ion_priority(ion) < _ion_priority(existing):
                self.peak_annotations[ion.exp_mz] = ion

        # ---------------------------------------------------------------
        # Step 6: Resolved matches for scoring (Phase 0 regression
        # hardening). Runs AFTER the Phase 0b envelope filter so scored
        # ions do not include precursor-envelope overlaps.
        # ---------------------------------------------------------------
        # `matched_ions` allows multiple theoretical ions to land on the
        # same experimental peak (e.g. a b-ion and y-ion that coincide,
        # or charge-reduced precursor envelope peaks that get multiple
        # glycan-containing c/z/Y labels). For coverage/FMR scoring we
        # need one ion per peak. Pick the closest-ppm winner.
        #
        # `matched_ions` is left unchanged for backward compatibility
        # with existing renderers that rely on full per-ion lists.
        scoring_resolved = {}
        for ion in self.matched_ions:
            prev = scoring_resolved.get(ion.exp_mz)
            if prev is None or abs(ion.mass_error_ppm) < abs(prev.mass_error_ppm):
                scoring_resolved[ion.exp_mz] = ion
        self.resolved_matches_for_scoring = list(scoring_resolved.values())

        # Build isotope-sibling map against the RAW spectrum (not deisotoped)
        # so envelopes that were collapsed during deisotoping are still
        # colored in the plot. For each matched ion, find raw-spectrum peaks
        # at mono ± k × (1.003355/charge) for k=1..4 and record them as
        # "isotope of" the parent. The render loop uses this to color those
        # peaks with the parent's ion color (instead of unassigned gray);
        # labels are suppressed to avoid clutter — the visible result is
        # that a full isotope envelope of Y(intact) / Y-ladder / c-ions
        # with large glycan reads as one colored ion.
        ISOTOPE_SPACING = 1.003355
        ISOTOPE_TOL_PPM = 10.0
        self.isotope_of = {}  # raw-spectrum isotope_mz -> parent MatchedIon
        raw_mz_arr = self.raw_mz
        for ion in self.matched_ions:
            if ion.charge <= 0:
                continue
            # Include Y(intact), charge-reduced precursor, c/z, b/y — any ion
            # with a detectable isotope envelope. Skipping only negative/zero
            # charge (which shouldn't happen but guards against bad data).
            mono = ion.exp_mz
            for k in (1, 2, 3, 4):
                target = mono + k * ISOTOPE_SPACING / ion.charge
                tol_da = target * ISOTOPE_TOL_PPM * 1e-6
                mask = np.abs(raw_mz_arr - target) <= tol_da
                if mask.any():
                    for idx in np.where(mask)[0]:
                        sib_mz = float(raw_mz_arr[idx])
                        # Don't override a peak that's already a real match
                        if sib_mz in self.matched_mz_set:
                            continue
                        # Keep the earliest-detected parent
                        if sib_mz not in self.isotope_of:
                            self.isotope_of[sib_mz] = ion

        # Find all glycan modification positions
        self.glycan_positions = self._detect_glycan_positions(modifications)

        # Calculate false match rate (spectrum shifting method from Schulte et al.)
        self.false_match_rate = calculate_false_match_rate(
            self.theoretical_ions, exp_mz, exp_intensity, tolerance_ppm
        )

        # Calculate comprehensive annotation statistics.
        # Pass the annotator's sn_threshold through so the stats function
        # uses the same filter the butterfly coverage view uses; otherwise
        # annotation_stats reports 0% coverage while the butterfly shows
        # real b/y matches (bug #10 in the reuse-failure writeup).
        self.annotation_stats = calculate_annotation_statistics(
            self.matched_ions, self.theoretical_ions,
            exp_mz, exp_intensity, len(peptide),
            min_sn=self.sn_threshold,
        )

    @staticmethod
    def _detect_glycan_positions(modifications: List[Dict]) -> List[int]:
        """Detect all glycan modification positions from the modifications list.

        Recognizes glycans by:
        1. Name containing glycan monosaccharide names
        2. Known glycan masses (528.2859, 203.0794, 299.123)
        3. Modifications on S/T with mass > 150 Da
        """
        glycan_name_parts = {'hex', 'hexnac', 'neuac', 'neugc', 'fuc',
                             'glyc', 'glcnac', 'galnac', 'sialyl'}
        known_glycan_masses = [528.2859, 203.0794, 299.123]
        positions = []

        for mod in modifications:
            pos = mod.get('position', 0)
            if pos <= 0:
                continue

            name_lower = mod.get('name', '').lower()
            if any(gn in name_lower for gn in glycan_name_parts):
                positions.append(pos)
                continue

            if any(abs(mod['mass'] - m) < 0.5 for m in known_glycan_masses):
                positions.append(pos)
                continue

            if mod.get('residue', '') in ('S', 'T') and mod['mass'] > 150:
                positions.append(pos)
                continue

        return sorted(positions)

    @staticmethod
    def _parse_glycan_label(label: str) -> Dict[str, int]:
        """Parse glycan label string into monosaccharide counts.

        Accepts both full and short formats:
        - Full: 'HexNAc(2)·Hex' -> {'HexNAc': 2, 'Hex': 1, ...}
        - Short: 'N2H1A2' or 'NHA' (implicit count=1) ->
                 {'HexNAc': 2, 'Hex': 1, 'NeuAc': 2, ...}

        Bare letters (no digits) are treated as count=1. Without this, multi-
        site glycopeptides where one site has a single-residue stub (e.g. 'N'
        for a lone HexNAc) silently drop that site from the combined-label
        aggregation, producing wrong per-fragment composition labels.
        """
        counts = {'HexNAc': 0, 'Hex': 0, 'NeuAc': 0, 'Fuc': 0}
        if not label:
            return counts
        # Full format FIRST when the label carries full monosaccharide names,
        # e.g. Byonic 'HexNAc(1)Hex(1)NeuAc(1)' (concatenated). Must precede the
        # short [NHFA] scan, else the H/N/A letters inside 'HexNAc'/'Hex'/'NeuAc'
        # get mis-read as short codes and a lone HexNAc(1) parses as N1H1A1.
        if re.search(r'HexNAc|NeuGc|NeuAc|dHex|Hex|Fuc', label):
            alias = {'HexNAc': 'HexNAc', 'Hex': 'Hex', 'Fuc': 'Fuc',
                     'dHex': 'Fuc', 'NeuAc': 'NeuAc', 'NeuGc': 'NeuAc'}
            any_match = False
            for m in re.finditer(r'(HexNAc|NeuGc|NeuAc|dHex|Hex|Fuc)(?:\((\d+)\))?', label):
                name = alias.get(m.group(1))
                if name:
                    counts[name] += int(m.group(2)) if m.group(2) else 1
                    any_match = True
            if any_match:
                return counts
        short_map = {'N': 'HexNAc', 'H': 'Hex', 'F': 'Fuc', 'A': 'NeuAc'}
        for m in re.finditer(r'([NHFA])(\d*)', label):
            name = short_map.get(m.group(1))
            if name:
                counts[name] += int(m.group(2)) if m.group(2) else 1
        return counts

    @staticmethod
    def _format_glycan_counts(counts: Dict[str, int]) -> str:
        """Format monosaccharide counts dict into short label string.

        e.g., {'HexNAc': 2, 'Hex': 1, ...} -> 'N2H1'
        """
        short = {'HexNAc': 'N', 'Hex': 'H', 'NeuAc': 'A', 'Fuc': 'F'}
        parts = []
        for name in ['HexNAc', 'Hex', 'Fuc', 'NeuAc']:
            n = counts.get(name, 0)
            if n > 0:
                parts.append(f"{short[name]}{n}")
        return ''.join(parts) if parts else "Glyc"

    def _combine_glycan_labels(self, positions: List[int]) -> str:
        """Combine glycan labels from multiple positions into one composition.

        Sums monosaccharide counts across all covered glycan sites.
        e.g., positions covering HexNAc·Hex + HexNAc -> HexNAc(2)·Hex
        """
        total = {'HexNAc': 0, 'Hex': 0, 'NeuAc': 0, 'Fuc': 0}
        for pos in positions:
            label = self.glycan_labels.get(pos, '')
            if label:
                counts = self._parse_glycan_label(label)
                for key in total:
                    total[key] += counts[key]
        return self._format_glycan_counts(total)

    def _get_ion_glycan_label(self, ion_type: str, ion_number: int) -> str:
        """Get glycan label string for a specific fragment ion.

        Determines which glycan site(s) the fragment covers and returns
        the combined composition label.

        Args:
            ion_type: 'b', 'c', 'y', or 'z'
            ion_number: Ion number (e.g., 6 for b6)

        Returns:
            Label string like '+HexNAc', '+HexNAc(2)·Hex', or ''
        """
        if not self.glycan_positions:
            return ""

        # Determine which glycan positions this ion covers
        covered = []
        if ion_type in ('b', 'c'):
            covered = [p for p in self.glycan_positions if p <= ion_number]
        elif ion_type in ('y', 'z'):
            covered = [p for p in self.glycan_positions
                       if p >= (len(self.peptide) - ion_number + 1)]

        if not covered:
            return ""

        if self.glycan_labels:
            combined = self._combine_glycan_labels(covered)
            return f"+{combined}"

        # No explicit labels — derive from modification masses at covered
        # positions. Fragment m/z is computed with the FULL glycan baked in
        # at each glycosite, so the label must reflect the same composition.
        # Falling back to "+HexNAc" was misleading when complex glycans were
        # actually attached (e.g., HexNAc2Hex1Fuc2 showed as "+HexNAc").
        from .fragment_calculator import FragmentCalculator
        total_mass = 0.0
        for mod in self.modifications:
            if mod.get('position', 0) in covered:
                total_mass += float(mod.get('mass', 0))
        if total_mass <= 0:
            return "+HexNAc"
        # OPair/FragPipe may encode a HexNAc-bearing O-glycan site as
        # 299.123 Da instead of the bare 203.0794 Da or TMT-tagged 528.2859 Da
        # mass used in other exports. These are representation choices for the
        # site modification, not distinct monosaccharide compositions, so keep
        # the reader-facing fragment label biological.
        if abs(total_mass - 299.123) < 0.5 or abs(total_mass - 203.0794) < 0.5:
            return "+HexNAc"
        if abs(total_mass - 528.2859) < 0.5:
            return "+HexNAc-TMT"

        decomp = FragmentCalculator._decompose_glycan_mass(total_mass, tol=0.05)
        if decomp is None:
            return f"+{total_mass:.2f}"
        nN, nH, nA, nF = decomp
        short = []
        if nN: short.append(f"N{nN}" if nN > 1 else "N")
        if nH: short.append(f"H{nH}" if nH > 1 else "H")
        if nA: short.append(f"A{nA}" if nA > 1 else "A")
        if nF: short.append(f"F{nF}" if nF > 1 else "F")
        return "+" + ("".join(short) if short else "Glyc")

    def _get_ion_color(self, ion: 'MatchedIon', has_neutral_loss: bool = False) -> str:
        """Get color for an ion type."""
        ion_type = ion.ion_type

        # Special handling for Y ions based on ion_number (Y0 vs Y1)
        if ion_type == 'Y':
            if ion.ion_number == 0:  # Y0 - glycan loss
                return ION_COLORS.get('Y0', ION_COLORS['Y'])
            else:  # Y1 - intact glycopeptide
                return ION_COLORS.get('Y1', ION_COLORS['Y'])

        # Site-determining ions (with glycan) always use solid colors
        if ion.has_modification:
            return ION_COLORS.get(ion_type, ION_COLORS['unassigned'])

        if has_neutral_loss:
            return ION_COLORS_NL.get(f'{ion_type}_NL', ION_COLORS.get(ion_type, ION_COLORS['unassigned']))
        return ION_COLORS.get(ion_type, ION_COLORS['unassigned'])

    def _format_annotation(self, ion: MatchedIon, short: bool = False) -> str:
        """Format ion annotation for display."""
        if ion.ion_type == 'oxonium':
            return ion.annotation
        elif ion.ion_type == 'precursor':
            # Format as [M+nH]+n for precursor ions
            if 'charge_reduced' in ion.annotation.lower() or 'cr' in ion.annotation.lower():
                return f"[M+{ion.charge}H]+{ion.charge}" if short else ion.annotation
            return f"[M+{ion.charge}H]+{ion.charge}" if short else ion.annotation
        elif ion.ion_type == 'Y':
            # Use the full composition label (e.g., Y0, Y0+N, Y0+NH, Y0+NHA)
            return ion.annotation
        else:
            # b, y, c, z ions
            charge_str = "" if ion.charge == 1 else f" {ion.charge}+"  # Add space before charge
            # Add glycan marker if ion contains glycan
            if ion.has_modification:
                glycan_str = self._get_ion_glycan_label(ion.ion_type, ion.ion_number)
            else:
                glycan_str = ""
            nl_str = ""
            if ion.neutral_loss:
                nl_map = {
                    'H2O': '°',
                    'NH3': '*',
                    'CO2': '**',
                    'HexNAc_TMT': '-Glyc',
                    'HexNAc': '-Hex',
                    # Bare b/y from complete glycan loss. Show the
                    # MSFragger-style backbone label (b5/y7) and keep the
                    # glycan-loss chemistry in the evidence table, not as a
                    # water-loss-like symbol.
                    'glyc_full': '',
                }
                nl_str = nl_map.get(ion.neutral_loss, f'-{ion.neutral_loss}')
            return f"{ion.ion_type}{ion.ion_number}{glycan_str}{nl_str}{charge_str}"

    @staticmethod
    def _ion_display_suppressed(ion) -> bool:
        """Whether an isotope-flagged ion should be hidden from the plot.

        MS2 spectra in this workflow are centroided peak lists. A missing or
        weak M+1 centroid is useful audit metadata, but it is too aggressive as
        a hard display filter: real glycan-retaining c/z ions and Y(intact)
        can be present as only the monoisotopic centroid after peak picking.

        Isotope flags are therefore used by peak-annotation priority to let an
        unflagged competing ion win on the same peak, but are not used to hide
        otherwise valid low-ppm matches outright.
        """
        return False

    def _get_fragmentation_coverage(self) -> Dict[str, set]:
        """Get which peptide bonds were fragmented by each ion type.

        Includes both base ions and neutral loss ions, since neutral loss
        ions also indicate bond cleavage.

        Uses S/N threshold (set in __init__) to prevent noise peaks from
        inflating coverage statistics. S/N is calculated from instrument
        noise if available, otherwise from median peak intensity.

        Phase 0h: ions that were display-suppressed by the isotope-
        consistency flag are also excluded here, so the sequence-arrow
        diagram and the "Fragmented Bonds X/Y" info line stay consistent
        with what is actually visible in the spectrum panel below.
        """
        coverage = {
            'b': set(),  # b ions (N-terminal, HCD)
            'c': set(),  # c ions (N-terminal, ETD)
            'y': set(),  # y ions (C-terminal, HCD)
            'z': set(),  # z ions (C-terminal, ETD)
        }

        for ion in self.matched_ions:
            if ion.ion_type not in coverage:
                continue
            if self._ion_display_suppressed(ion):
                continue
            ion_sn = getattr(ion, '_sn', float('inf'))
            if ion_sn >= self.sn_threshold:
                coverage[ion.ion_type].add(ion.ion_number)

        return coverage

    def _get_glycan_containing_ions(self) -> Dict[str, set]:
        """Get which matched ions contain the glycan (site-determining ions).

        Returns a dict mapping ion type to set of ion numbers that contain glycan.

        Phase 0h: display-suppressed ions are excluded, matching the
        spectrum display.
        """
        glycan_ions = {
            'b': set(),
            'c': set(),
            'y': set(),
            'z': set(),
        }

        for ion in self.matched_ions:
            if ion.ion_type not in glycan_ions:
                continue
            if self._ion_display_suppressed(ion):
                continue
            if ion.has_modification:
                glycan_ions[ion.ion_type].add(ion.ion_number)

        return glycan_ions

    def plot_peptide_butterfly(self,
                               ax=None,
                               figsize: Tuple[float, float] = (4.0, 1.2),
                               fontsize: Optional[dict] = None,
                               output_path: Optional[str] = None,
                               coverage: Optional[Dict[str, set]] = None):
        """The IPSA-style peptide coverage ladder — b/c marks above, y/z below — on its own.

        EXTRACTED FROM ``plot()`` 2026-08-23, unchanged. ``plot()`` now calls this for its section
        1, so there is one implementation and the composed figure cannot drift from the standalone
        panel. It was ~170 lines inline and reachable only by rendering the whole 3-4 row figure.

        WHY IT IS PUBLIC AND TAKES AN ``ax``. Callers wanting only the coverage ladder — a figure
        panel, a review sheet, a slide — had to call ``plot()`` and discard three quarters of what
        it returned, inheriting its GridSpec, its font scheme and its ``tight_layout()``. Those are
        reasonable defaults for a review PDF and wrong for a manuscript panel with its own
        standard. Passing ``ax`` draws into someone else's figure and returns it untouched
        otherwise.

        Args:
            ax: draw into this Axes. When None a new figure of ``figsize`` is created and returned.
            figsize: only used when ``ax`` is None.
            fontsize: override the two sizes this panel uses, ``{'seq': pt, 'frag': pt}``. Defaults
                to ``plot(compact=True)``'s values (7 and 5). A caller with a font floor above 5 pt
                — the Nature-family 7 pt minimum, for instance — needs this; there was no way to
                reach it before.
            output_path: save the figure as PDF. Ignored when ``ax`` was supplied, since the
                figure then belongs to the caller.
            coverage: override the covered-bond sets, ``{'b': {..}, 'c': set(), 'y': {..},
                'z': set()}``. By default the ladder draws this annotator's OWN matching
                (``_get_fragmentation_coverage()``: S/N-gated, neutral losses count). A caller
                placing the ladder under a spectrum panel annotated by a DIFFERENT pipeline
                passes that pipeline's covered bonds here, so the two stacked panels cannot
                disagree about which ions exist. Missing keys are treated as empty.

        Returns:
            ``(fig, coverage)`` — ``fig`` is None when ``ax`` was supplied. ``coverage`` is the
            per-ion-type covered-position dict, which the caller needs for a coverage statistic and
            which ``plot()`` goes on to use for its info panel.
        """
        FZ = dict(seq=7, frag=5)
        if fontsize:
            FZ.update(fontsize)

        fig = None
        if ax is None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        ax_seq = ax

        # =====================================================================
        # 1. Peptide Sequence with Fragmentation Sites (IPSA-style)
        # =====================================================================
        ax_seq.set_xlim(-0.5, len(self.peptide) + 0.5)
        # ylim will be adjusted after labels are placed
        ax_seq.axis('off')

        # Get fragmentation coverage by ion type -- the caller's set when supplied, else ours
        if coverage is None:
            coverage = self._get_fragmentation_coverage()
        else:
            coverage = {t: set(coverage.get(t) or ()) for t in ('b', 'c', 'y', 'z')}
        # Get glycan-containing (site-determining) ions
        glycan_ions = self._get_glycan_containing_ions()

        # Draw peptide sequence (no boxes, Arial font, size 8)
        # Collect arrow labels for stagger placement
        _arrow_labels_top = []   # (x, base_y, text, color) for N-terminal (b,c)
        _arrow_labels_bot = []   # (x, base_y, text, color) for C-terminal (y,z)

        for i, aa in enumerate(self.peptide):
            pos = i + 1  # 1-based
            # Highlight glycosylated residues with color only (no box)
            if pos in self.glycan_positions:
                text_color = '#D2691E'  # Chocolate color for glycosylated residue
                fontweight = 'bold'
            elif pos in self.n_sequon_positions or (pos - 1) in self.n_sequon_positions or (pos - 2) in self.n_sequon_positions:
                # N-X-S/T sequon: highlight N and the two following residues in purple
                text_color = '#8B008B'  # Dark magenta for N-sequon
                fontweight = 'bold'
            else:
                text_color = 'black'
                fontweight = 'normal'

            ax_seq.text(i + 0.5, 0.5, aa, ha='center', va='center', fontsize=FZ['seq'],
                       fontweight=fontweight, fontfamily='Arial', color=text_color)

            # Glycan composition label below residue — disabled for cleaner plots
            # if (i + 1) in self.glycan_labels:
            #     gl = self.glycan_labels[i + 1]
            #     ax_seq.text(i + 0.5, -0.15, gl, ...)

            # Draw IPSA-style fragmentation marks (vertical + diagonal)
            # b/y have longer vertical, c/z have shorter vertical, same x position
            if i < len(self.peptide) - 1:
                bond_pos = i + 1  # N-terminal ion number (b1, c1 after first AA)
                c_bond = len(self.peptide) - bond_pos  # C-terminal ion number

                # Base position for marks (between amino acids)
                x_base = i + 0.85

                # N-terminal ions (b/c) - vertical from middle up, then diagonal up-LEFT
                # b ions (blue) - longer vertical
                if bond_pos in coverage['b']:
                    y_top_b = 0.9
                    lw = 2.5 if bond_pos in glycan_ions['b'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_top_b], color=ION_COLORS['b'], linewidth=lw)
                    ax_seq.plot([x_base, x_base - 0.15], [y_top_b, y_top_b + 0.10], color=ION_COLORS['b'], linewidth=lw)
                    # Every covered bond gets its ion label; the glycan tag joins it only when
                    # the ion carries one. (Labels were glycan-only before, which left a plain
                    # peptide's ladder as bare ticks with nothing naming them.)
                    g_tag = self._get_ion_glycan_label('b', bond_pos) if bond_pos in glycan_ions['b'] else ''
                    _arrow_labels_top.append((x_base - 0.6, y_top_b + 0.25,
                                              f'←b{bond_pos}{g_tag}', ION_COLORS['b']))
                # c ions (green) - shorter vertical
                if bond_pos in coverage['c']:
                    y_top_c = 0.75
                    lw = 2.5 if bond_pos in glycan_ions['c'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_top_c], color=ION_COLORS['c'], linewidth=lw)
                    ax_seq.plot([x_base, x_base - 0.15], [y_top_c, y_top_c + 0.10], color=ION_COLORS['c'], linewidth=lw)
                    g_tag = self._get_ion_glycan_label('c', bond_pos) if bond_pos in glycan_ions['c'] else ''
                    _arrow_labels_top.append((x_base - 0.6, y_top_c + 0.15,
                                              f'←c{bond_pos}{g_tag}', ION_COLORS['c']))

                # C-terminal ions (y/z) - vertical from middle down, then diagonal down-RIGHT
                # y ions (red) - longer vertical
                if c_bond in coverage['y']:
                    y_bot_y = 0.1
                    lw = 2.5 if c_bond in glycan_ions['y'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_bot_y], color=ION_COLORS['y'], linewidth=lw)
                    ax_seq.plot([x_base, x_base + 0.15], [y_bot_y, y_bot_y - 0.10], color=ION_COLORS['y'], linewidth=lw)
                    g_tag = self._get_ion_glycan_label('y', c_bond) if c_bond in glycan_ions['y'] else ''
                    _arrow_labels_bot.append((x_base + 0.2, y_bot_y - 0.30,
                                              f'y{c_bond}{g_tag}→', ION_COLORS['y']))
                # z ions (orange) - shorter vertical
                if c_bond in coverage['z']:
                    y_bot_z = 0.25
                    lw = 2.5 if c_bond in glycan_ions['z'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_bot_z], color=ION_COLORS['z'], linewidth=lw)
                    ax_seq.plot([x_base, x_base + 0.15], [y_bot_z, y_bot_z - 0.10], color=ION_COLORS['z'], linewidth=lw)
                    g_tag = self._get_ion_glycan_label('z', c_bond) if c_bond in glycan_ions['z'] else ''
                    _arrow_labels_bot.append((x_base + 0.2, y_bot_z - 0.20,
                                              f'z{c_bond}{g_tag}→', ION_COLORS['z']))

        # Deduplicate arrow labels: when consecutive ions of the SAME TYPE carry
        # the same glycan tag, keep only the first label (arrows still drawn).
        # Different ion types (b vs c, y vs z) always get their own labels.
        def _dedup_labels(labels_list):
            if not labels_list:
                return labels_list
            labels_list.sort(key=lambda l: l[0])
            kept = [labels_list[0]]
            for entry in labels_list[1:]:
                text = entry[2]
                prev_text = kept[-1][2]
                # Extract ion type letter (b, c, y, z) from label text
                import re as _re
                cur_type = _re.search(r'[bcyz]', text)
                prev_type = _re.search(r'[bcyz]', prev_text)
                cur_type = cur_type.group() if cur_type else ''
                prev_type = prev_type.group() if prev_type else ''
                # Only dedup if same ion type
                if cur_type != prev_type:
                    kept.append(entry)
                    continue
                # Extract glycan tag (part after "+")
                cur_tag = text.split('+', 1)[1] if '+' in text else ''
                prev_tag = prev_text.split('+', 1)[1] if '+' in prev_text else ''
                # Strip arrow symbols for comparison
                cur_tag = cur_tag.rstrip('→').lstrip('←')
                prev_tag = prev_tag.rstrip('→').lstrip('←')
                if cur_tag and cur_tag == prev_tag:
                    continue  # skip duplicate glycan label of same ion type
                kept.append(entry)
            return kept

        _arrow_labels_top = _dedup_labels(_arrow_labels_top)
        _arrow_labels_bot = _dedup_labels(_arrow_labels_bot)

        # Place arrow labels with stagger to avoid overlap.
        #
        # The extents are estimated in PHYSICAL units from the figure's real geometry and then
        # converted to data coordinates. The previous constant (0.06 * (len+1) / 14) was a
        # data-coordinate guess calibrated on the default 4.0 in figure: on a narrower figure it
        # under-estimated label widths by the width ratio, the stagger never fired, and a dense
        # ladder (every covered bond labelled, 2026-08-23) rendered its labels on top of each
        # other. 0.62 em per character is Arial's average advance; 1.3 em is the line height.
        _fig_ref = ax_seq.figure
        _bbox = ax_seq.get_position()
        _axes_w_in = max(_fig_ref.get_size_inches()[0] * _bbox.width, 1e-6)
        _axes_h_in = max(_fig_ref.get_size_inches()[1] * _bbox.height, 1e-6)
        _x_range = len(self.peptide) + 1.0
        char_width = (FZ['frag'] * 0.62 / 72.0) / _axes_w_in * _x_range
        # Vertical step between staggered rows, sized to the label's line height against the
        # ladder's ~2.0-unit baseline y-range (the final ylim grows with the stack, which only
        # makes the physical step conservative).
        stagger_step = max(0.18, (FZ['frag'] * 1.3 / 72.0) / _axes_h_in * 2.0)

        def _stagger_rows(labels_list, direction):
            """Greedy ROW assignment: a label moves to the next row only while the row it is
            trying holds an overlapping label. Comparing against every placed label regardless
            of row (the previous behaviour) turned a dense ladder into an ever-deepening
            staircase -- each label fleeing the one before it -- when two alternating rows hold
            the same set comfortably."""
            placed = []  # (x, y, est_width)
            for x, base_y, text, color in sorted(labels_list, key=lambda l: l[0]):
                est_width = len(text) * char_width
                y = base_y
                while any(abs(x - px) < max(est_width, pw) * 0.85
                          and abs(y - py) < stagger_step * 0.5
                          for px, py, pw in placed):
                    y = y + stagger_step if direction == 'up' else y - stagger_step
                placed.append((x, y, est_width))
                yield x, y, text, color

        for labels_list, direction in [(_arrow_labels_top, 'up'), (_arrow_labels_bot, 'down')]:
            for x, y, text, color in _stagger_rows(labels_list, direction):
                ax_seq.text(x, y, text, fontsize=FZ['frag'], color=color,
                           fontweight='bold', fontfamily='Arial')

        # Set ylim dynamically based on placed labels and glycan labels
        y_top = 1.5
        y_bot = -0.5
        if _arrow_labels_top:
            ys = [y for _x, y, _t, _c in _stagger_rows(_arrow_labels_top, 'up')]
            if ys:
                y_top = max(y_top, max(ys) + 0.3)
        if _arrow_labels_bot:
            # Bottom labels exist on every covered y/z bond now, and a staggered stack must not
            # fall off the axis.
            ys = [y for _x, y, _t, _c in _stagger_rows(_arrow_labels_bot, 'down')]
            if ys:
                y_bot = min(y_bot, min(ys) - 0.3)
        ax_seq.set_ylim(y_bot, y_top)

        if fig is not None and output_path:
            fig.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
        return fig, coverage

    def plot(self,
             figsize: Tuple[float, float] = (6, 4),
             output_path: Optional[str] = None,
             show_error_plot: bool = True,
             intensity_threshold_pct: float = 1.0,
             max_labels: int = 50,
             compact: bool = True,
             mini: bool = False) -> plt.Figure:
        """
        Create annotated spectrum plot.

        Args:
            figsize: Figure size (width, height)
            output_path: Path to save PDF (optional)
            show_error_plot: Whether to show mass error plot
            intensity_threshold_pct: Minimum intensity (% base peak) to label
            max_labels: Maximum number of labels to show

        Returns:
            matplotlib Figure object
        """
        # Set Arial as default font for all text
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.sans-serif'] = ['Arial']
        # Font-size scheme. compact=True targets small publication panels
        # (~6x4 in, 7-8 pt main text; dense peak/fragment labels a touch
        # smaller so they still fit). Default preserves the 12x8 hierarchy.
        if mini:
            FZ = dict(seq=5, frag=4, info1=6, info2=5, peak=5, thresh=4,
                      axis=6, tick=5, legend=5, symleg=5, title=[7, 6, 6, 6],
                      suptitle=7, subfile=6)
        elif compact:
            FZ = dict(seq=7, frag=5, info1=8, info2=7, peak=6, thresh=5,
                      axis=8, tick=7, legend=6, symleg=6, title=[8, 7, 7, 7],
                      suptitle=8, subfile=7)
        else:
            FZ = dict(seq=8, frag=5, info1=9, info2=8, peak=7, thresh=6,
                      axis=9, tick=8, legend=7, symleg=6, title=[12, 10, 10, 10],
                      suptitle=10, subfile=8)

        # Normalize intensities to % base peak.
        # IMPORTANT: the PLOT uses the RAW spectrum so full isotope envelopes
        # of matched ions remain visible — even ones that aggressive
        # deisotoping collapsed. Matching itself was done on the deisotoped
        # self.exp_mz/self.exp_intensity to avoid false matches at isotope
        # m/z positions. The base peak here is the raw base peak so relative
        # intensities read normally.
        base_peak = float(np.max(self.raw_intensity))
        plot_mz = self.raw_mz
        plot_intensity = self.raw_intensity
        rel_intensity = plot_intensity / base_peak * 100

        # Set up figure
        seq_height = 1.2
        if mini:
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(2, 1, height_ratios=[seq_height, 3.4], hspace=0.04)
            ax_seq = fig.add_subplot(gs[0])
            ax_spec = fig.add_subplot(gs[1])
            ax_info = None
            ax_error = None
        elif show_error_plot:
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(4, 1, height_ratios=[seq_height, 0.5, 3, 1], hspace=0.05)
            ax_seq = fig.add_subplot(gs[0])
            ax_info = fig.add_subplot(gs[1])
            ax_spec = fig.add_subplot(gs[2])
            ax_error = fig.add_subplot(gs[3], sharex=ax_spec)
        else:
            fig = plt.figure(figsize=(figsize[0], figsize[1] * 0.8))
            gs = GridSpec(3, 1, height_ratios=[seq_height, 0.5, 3], hspace=0.05)
            ax_seq = fig.add_subplot(gs[0])
            ax_info = fig.add_subplot(gs[1])
            ax_spec = fig.add_subplot(gs[2])
            ax_error = None

        # =====================================================================
        # 1. Peptide Sequence with Fragmentation Sites (IPSA-style)
        # =====================================================================
        # Extracted to plot_peptide_butterfly() 2026-08-23 so the ladder is reachable on its own.
        # `coverage` is returned because the info panel below needs it.
        _, coverage = self.plot_peptide_butterfly(
            ax=ax_seq, fontsize={'seq': FZ['seq'], 'frag': FZ['frag']})

        # =====================================================================
        # 2. Info Panel (with False Match Rate)
        # =====================================================================
        if ax_info is not None:
            ax_info.axis('off')

        # Calculate statistics - convert ion numbers to bond positions
        total_bonds = len(self.peptide) - 1
        covered_bonds = set()
        for pos in (coverage['b'] | coverage['c']):
            covered_bonds.add(pos)  # b_n/c_n covers bond n
        for pos in (coverage['y'] | coverage['z']):
            covered_bonds.add(len(self.peptide) - pos)  # y_n/z_n covers bond (len-n)
        fragmented_bonds = len(covered_bonds)
        coverage_pct = fragmented_bonds / total_bonds * 100 if total_bonds > 0 else 0

        # Get false match rate values
        fmr = self.false_match_rate
        stats = self.annotation_stats

        # Line 1: Basic info with fragmented bonds (matches IPSA reference style)
        info_line1 = (f"Gene: {self.gene}    Site: {self.site_index}    "
                     f"Precursor m/z: {self.precursor_mz:.4f}    "
                     f"Charge: +{self.precursor_charge}    "
                     f"Fragmented Bonds: {fragmented_bonds}/{total_bonds} ({coverage_pct:.0f}%)")

        # Line 2: FMR
        info_line2 = f"FMR: {fmr.fmr_peaks*100:.1f}%"

        if ax_info is not None:
            ax_info.text(0.5, 0.75, info_line1, ha='center', va='center', fontsize=FZ['info1'],
                        fontfamily='Arial', transform=ax_info.transAxes)
            ax_info.text(0.5, 0.25, info_line2, ha='center', va='center', fontsize=FZ['info2'],
                        fontfamily='Arial', transform=ax_info.transAxes, color='#555555')

        # =====================================================================
        # 3. Spectrum Plot
        # =====================================================================
        # Phase 0g: totally suppress display for isotope-flag-failed ions
        # (M1_absent, suspicious_low_M1). Previously Phase 0f only dropped
        # the LABEL; the colored vline + error-plot dot still drew, making
        # the peak appear matched (wrong color) with no label. Now the peak
        # reverts to gray/unmatched appearance and the error-plot dot goes
        # away too. Reviewer rule still holds: ion stays in matched_ions and
        # resolved_matches_for_scoring — this is DISPLAY ONLY.
        _is_display_suppressed = self._ion_display_suppressed

        # Peaks whose ONLY matched ion(s) are all display-suppressed must
        # render as gray unmatched peaks (instead of being invisible).
        from collections import defaultdict as _dd
        _ions_by_mz = _dd(list)
        for _ion in self.matched_ions:
            _ions_by_mz[_ion.exp_mz].append(_ion)
        _suppressed_only_mz = {
            _mz for _mz, _ions in _ions_by_mz.items()
            if _ions and all(_is_display_suppressed(_i) for _i in _ions)
        }

        # Plot all peaks first (from RAW spectrum so isotope envelopes are
        # visible). Unmatched peaks get the unassigned color, except those
        # identified as isotope siblings of a matched ion — those inherit
        # the parent's color (no separate label) so full isotope envelopes
        # of Y(intact) / Y-ladder / c-ions-with-glycan read as one ion.
        for mz, inten in zip(plot_mz, rel_intensity):
            # Treat suppressed-only peaks as unmatched (fall through to gray)
            if mz in self.matched_mz_set and mz not in _suppressed_only_mz:
                continue  # drawn in the matched loop below with proper color
            parent = self.isotope_of.get(mz)
            if parent is not None and not _is_display_suppressed(parent):
                iso_color = self._get_ion_color(
                    parent,
                    bool(parent.neutral_loss and parent.neutral_loss != 'glyc_full'),
                )
                ax_spec.vlines(mz, 0, inten, color=iso_color, linewidth=1.0, alpha=0.75)
            else:
                ax_spec.vlines(mz, 0, inten, color=ION_COLORS['unassigned'], linewidth=0.8)

        # ---------------------------------------------------------------
        # Phase 0d fix: draw colored vlines from ALL matched_ions (unchanged),
        # but place LABELS only for the peak_annotations winner per peak.
        # The old label loop used its own `label_priority`
        # (has_modification, -intensity) which DID NOT respect isotope_flag
        # or the Patch D `_ion_priority`. That caused Y0 2+ to beat y14 1+ on
        # the same peak in case 5: both had has_modification=False AND the
        # exact same exp_intensity (same peak!), so the tuple tied and Y0's
        # original insertion order won. peak_annotations was already resolving
        # to y14 1+ (per Patch D), but the plot never used it. Now it does.
        # ---------------------------------------------------------------
        # Step A: draw colored peak lines for all matched ions. The same
        # exp_mz may be painted multiple times (one per matched ion); the last
        # draw wins visibility. We keep the original sort so glycan-bearing
        # ions are drawn first and non-glycan ions last — i.e. on a peak with
        # a Y0 and a backbone y ion, the visible color matches the backbone
        # winner y14 (red), which is consistent with the y14 label we will
        # place below.
        def _vline_order(ion):
            priority = 0 if ion.has_modification else 1
            return (priority, -ion.exp_intensity)
        for ion in sorted(self.matched_ions, key=_vline_order):
            if ion.ion_type == 'precursor':
                continue
            if _is_display_suppressed(ion):
                continue  # Phase 0g: skip colored vline for suppressed ions
            color = self._get_ion_color(
                ion,
                bool(ion.neutral_loss and ion.neutral_loss != 'glyc_full'),
            )
            rel_int = ion.exp_intensity / base_peak * 100
            ax_spec.vlines(ion.exp_mz, 0, rel_int, color=color, linewidth=1.5)

        # Step B: label only one ion per peak, using the Patch-D-resolved
        # winner in self.peak_annotations. This is where isotope_flag takes
        # effect on what the reader actually sees.
        # The geometric policy -- intensity order, budget, floor, the 15 m/z x 8 % exclusion,
        # rotation -- lives in resolve_peak_labels() (extracted 2026-08-23, pure). This loop
        # keeps only the ION-specific part: precursor/isotope filtering, colour, formatting.
        cand = []
        style = {}
        for ion in self.peak_annotations.values():
            if ion.ion_type == 'precursor':
                continue
            # Isotope flags are audit/priority metadata for centroided MS2,
            # not a hard label filter. A missing M+1 centroid can occur after
            # peak picking even for real glycan-retaining c/z ions or
            # Y(intact). Only hide labels when an observed M+1 partner is
            # abnormally weak relative to the mono peak.
            flag = getattr(ion, 'isotope_flag', None)
            if flag == 'suspicious_low_M1':
                continue
            rel_int = ion.exp_intensity / base_peak * 100
            label = self._format_annotation(ion, short=True)
            cand.append((ion.exp_mz, rel_int, label))
            style[(float(ion.exp_mz), label)] = self._get_ion_color(
                ion,
                bool(ion.neutral_loss and ion.neutral_loss != 'glyc_full'),
            )
        for mz, rel_int, label, rot in resolve_peak_labels(
                cand, max_labels=max_labels,
                intensity_threshold_pct=intensity_threshold_pct):
            ax_spec.annotate(label, (mz, rel_int),
                            textcoords="offset points", xytext=(0, 5),
                            ha='center', va='bottom', fontsize=FZ['peak'],
                            fontfamily='Arial', color=style[(mz, label)], fontweight='bold',
                            rotation=rot)

        # 5% line only for HCD (trigger threshold), not for EThcD/ETD
        if self.activation_type == "HCD":
            ax_spec.axhline(5.0, color='gray', linewidth=0.8, linestyle='--', zorder=0)
            ax_spec.text(max(plot_mz) + 45, 5.0, '5%', fontsize=FZ['thresh'], fontfamily='Arial',
                        color='gray', va='center', ha='right')

        # Set spectrum axes
        ax_spec.set_xlim(min(plot_mz) - 50, max(plot_mz) + 50)
        ax_spec.set_ylim(0, 110)
        ax_spec.set_ylabel('Relative Abundance (%)', fontsize=FZ['axis'], fontfamily='Arial')
        ax_spec.spines['top'].set_visible(False)
        ax_spec.spines['right'].set_visible(False)
        ax_spec.tick_params(labelsize=FZ['tick'])

        if ax_error is None:  # no ppm panel below (show_error_plot off, or mini) -> label m/z here
            ax_spec.set_xlabel('m/z', fontsize=FZ['axis'], fontfamily='Arial')

        # =====================================================================
        # 4. Mass Error Plot
        # =====================================================================
        if ax_error is not None:
            for ion in self.matched_ions:
                if ion.ion_type == 'precursor':
                    continue
                if _is_display_suppressed(ion):
                    continue  # Phase 0g: no ppm-error scatter for suppressed ions
                color = self._get_ion_color(
                    ion,
                    bool(ion.neutral_loss and ion.neutral_loss != 'glyc_full'),
                )
                ax_error.scatter(ion.exp_mz, ion.mass_error_ppm, color=color, s=20, alpha=0.7)

            ax_error.axhline(0, color='black', linewidth=0.5, linestyle='--')
            ax_error.axhline(self.tolerance_ppm, color='red', linewidth=0.5, linestyle=':')
            ax_error.axhline(-self.tolerance_ppm, color='red', linewidth=0.5, linestyle=':')

            ax_error.set_ylim(-self.tolerance_ppm * 1.5, self.tolerance_ppm * 1.5)
            ax_error.set_xlabel('m/z', fontsize=FZ['axis'], fontfamily='Arial')
            ax_error.set_ylabel('Error (ppm)', fontsize=FZ['axis'], fontfamily='Arial')
            ax_error.spines['top'].set_visible(False)
            ax_error.spines['right'].set_visible(False)
            ax_error.tick_params(labelsize=FZ['tick'])

        # =====================================================================
        # 5. Legend
        # =====================================================================
        legend_elements = [
            Line2D([0], [0], color=ION_COLORS['b'], linewidth=2, label='b'),
            Line2D([0], [0], color=ION_COLORS['y'], linewidth=2, label='y'),
            Line2D([0], [0], color=ION_COLORS['c'], linewidth=2, label='c'),
            Line2D([0], [0], color=ION_COLORS['z'], linewidth=2, label='z'),
            Line2D([0], [0], color=ION_COLORS['Y0'], linewidth=2, label='Y0'),
            Line2D([0], [0], color=ION_COLORS['Y1'], linewidth=2, label='Y-ladder'),
            Line2D([0], [0], color=ION_COLORS['oxonium'], linewidth=2, label='Oxonium'),
        ]
        if not mini:
         ax_spec.legend(handles=legend_elements, loc='upper right', fontsize=FZ['legend'], ncol=4,
                      prop={'family': 'Arial', 'size': FZ['legend']}, handlelength=1.2, handletextpad=0.4,
                      columnspacing=0.8, borderpad=0.4)

        # Add annotation symbol legend (neutral losses and glycan)
        if self.glycan_labels:
            symbol_legend = "+N/H/F/A: with glycan (N=HexNAc H=Hex F=Fuc A=NeuAc)    °: -H2O    *: -NH3"
        else:
            symbol_legend = "+HexNAc: with glycan    °: -H2O    *: -NH3"
        if not mini:
         ax_spec.text(0.02, 0.97, symbol_legend, transform=ax_spec.transAxes,
                    fontsize=FZ['symleg'], fontfamily='Arial', color='#555555',
                    verticalalignment='top', horizontalalignment='left')

        # Title — caller may override with a composed multi-line title
        # (e.g. to include file | scan | activation | SA | glycan composition).
        # Use stacked fig.text placements to avoid matplotlib's suptitle
        # auto-wrap colliding with subsequent fig.text lines.
        if self.title_override:
            lines = self.title_override.split('\n')
            # Reserve ~15% of figure height for the 3-4 line title so the
            # spectrum axes never crash into the bottom title line.
            fig.subplots_adjust(top=0.90 if mini else 0.85)
            # Wide line spacing (0.035 of figure height ≈ 27 pt) ensures no
            # visual overlap even at 12 pt bold top line.
            y_positions = [0.965] if mini else [0.980, 0.945, 0.910, 0.875]
            font_sizes = [FZ['suptitle']] if mini else FZ['title']
            weights = ['bold', 'normal', 'normal', 'normal']
            colors = ['black', '#222222', '#333333', '#555555']
            for i, line in enumerate(lines[:len(y_positions)]):
                fig.text(0.5, y_positions[i], line, ha='center',
                         fontsize=font_sizes[i], fontweight=weights[i],
                         fontfamily='Arial', color=colors[i])
        else:
            title_parts = [self.gene, self.site_index]
            if self.scan_num is not None:
                title_parts.append(f's{self.scan_num}')
            if self.confidence_level:
                title_parts.append(self.confidence_level)
            act_display = {'ETHCD': 'EThcD', 'HCD': 'HCD', 'ETD': 'ETD', 'CID': 'CID'}
            title_parts.append(act_display.get(self.activation_type, self.activation_type))
            fig.suptitle(' - '.join(p for p in title_parts if p),
                        fontsize=FZ['suptitle'], fontweight='bold', fontfamily='Arial', y=0.98)
            if self.source_file:
                fig.text(0.5, 0.935, self.source_file, ha='center',
                         fontsize=FZ['subfile'], fontfamily='Arial', color='#666666')

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
            print(f"Saved: {output_path}")

        return fig


# =============================================================================
# BATCH ANNOTATION FUNCTION
# =============================================================================

def annotate_spectra_batch(summary_file: str,
                          spectra_dir: str,
                          output_dir: str,
                          n_spectra: int = None,
                          tolerance_ppm: float = 20.0,
                          save_statistics: bool = True) -> Tuple[List[str], Optional[pd.DataFrame]]:
    """
    Batch annotate multiple spectra.

    Args:
        summary_file: Path to extracted spectra summary CSV
        spectra_dir: Directory containing spectrum CSV files
        output_dir: Directory for output PDFs
        n_spectra: Number of spectra to process (None = all)
        tolerance_ppm: Mass tolerance in ppm
        save_statistics: Whether to save annotation statistics CSV

    Returns:
        Tuple of (list of output file paths, DataFrame with statistics)
    """
    # Load summary
    df = pd.read_csv(summary_file)

    if n_spectra:
        df = df.head(n_spectra)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    output_files = []
    statistics_list = []

    for idx, row in df.iterrows():
        print(f"\nProcessing {idx + 1}/{len(df)}: {row['site_index']}")

        try:
            # Load spectrum data
            spec_file = os.path.join(spectra_dir, row['spectrum_file'])
            spec_df = pd.read_csv(spec_file)
            exp_mz = spec_df['mz'].values
            exp_intensity = spec_df['intensity'].values

            # Parse modifications
            if 'modifications_json' in row and pd.notna(row['modifications_json']):
                modifications = json.loads(row['modifications_json'])
            else:
                modifications = parse_modifications_from_string(row.get('Assigned_Modifications', ''))

            # Create annotator
            annotator = SpectrumAnnotator(
                peptide=row['Peptide'],
                modifications=modifications,
                precursor_charge=int(row['Charge']),
                precursor_mz=float(row.get('Calibrated_Observed_MZ', row.get('mzML_precursor_mz', 0))),
                exp_mz=exp_mz,
                exp_intensity=exp_intensity,
                tolerance_ppm=tolerance_ppm,
                site_index=row['site_index'],
                gene=row['Gene']
            )

            # Generate output filename
            output_file = os.path.join(output_dir, f"{row['site_index']}_{row['scan_number']}.pdf")

            # Create plot
            fig = annotator.plot(output_path=output_file)
            plt.close(fig)

            output_files.append(output_file)

            # Collect statistics
            fmr = annotator.false_match_rate
            stats = annotator.annotation_stats
            statistics_list.append({
                'site_index': row['site_index'],
                'gene': row['Gene'],
                'peptide': row['Peptide'],
                'charge': int(row['Charge']),
                'scan_number': row.get('scan_number', ''),
                'sequence_coverage': stats['sequence_coverage'],
                'sequence_coverage_bonds': stats['sequence_coverage_bonds'],
                'peaks_annotated': stats['peaks_annotated'],
                'peaks_annotated_count': stats['peaks_annotated_count'],
                'intensity_annotated': stats['intensity_annotated'],
                'fragments_found': stats['fragments_found'],
                'fragments_found_count': stats['fragments_found_count'],
                'fmr_peaks': fmr.fmr_peaks,
                'fmr_intensity': fmr.fmr_intensity,
                'matched_peaks': fmr.matched_peaks,
                'matched_intensity': fmr.matched_intensity,
                'avg_random_peaks': fmr.avg_random_peaks,
                'avg_random_intensity': fmr.avg_random_intensity,
            })

        except Exception as e:
            print(f"  Error: {e}")
            continue

    # Create statistics DataFrame
    stats_df = pd.DataFrame(statistics_list) if statistics_list else None

    # Save statistics if requested
    if save_statistics and stats_df is not None:
        stats_file = os.path.join(output_dir, "annotation_statistics.csv")
        stats_df.to_csv(stats_file, index=False)
        print(f"\nSaved annotation statistics to: {stats_file}")

    return output_files, stats_df


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Test with 10 example spectra from HEK293T
    SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
    OUTPUT_PATH = os.path.join(SOURCE_PATH, "point_to_point_response/")

    summary_file = os.path.join(OUTPUT_PATH, "extracted_spectra_EThcD_HEK293T_summary.csv")
    spectra_dir = os.path.join(OUTPUT_PATH, "extracted_spectra_EThcD_HEK293T")
    output_dir = os.path.join(OUTPUT_PATH, "annotated_spectra_test")

    print("=" * 70)
    print("Spectrum Annotation Test - 10 Example Spectra")
    print("(with False Match Rate calculation)")
    print("=" * 70)

    output_files, stats_df = annotate_spectra_batch(
        summary_file=summary_file,
        spectra_dir=spectra_dir,
        output_dir=output_dir,
        n_spectra=10,
        tolerance_ppm=20.0,
        save_statistics=True
    )

    print("\n" + "=" * 70)
    print(f"Generated {len(output_files)} annotated spectra")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

    # Display summary statistics
    if stats_df is not None and len(stats_df) > 0:
        print("\nAnnotation Statistics Summary:")
        print("-" * 50)
        print(f"Mean Sequence Coverage: {stats_df['sequence_coverage'].mean()*100:.1f}%")
        print(f"Mean Intensity Annotated: {stats_df['intensity_annotated'].mean()*100:.1f}%")
        print(f"Mean FMR (peaks): {stats_df['fmr_peaks'].mean()*100:.1f}%")
        print(f"Mean FMR (intensity): {stats_df['fmr_intensity'].mean()*100:.1f}%")
        print("-" * 50)
        print("\nPer-spectrum statistics:")
        print(stats_df[['site_index', 'sequence_coverage', 'intensity_annotated',
                       'fmr_peaks', 'fmr_intensity']].to_string(index=False))
