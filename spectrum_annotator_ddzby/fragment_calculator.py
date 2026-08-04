#!/usr/bin/env python3
"""
Theoretical Fragment Ion Calculator for Glycopeptides

This module calculates theoretical m/z values for peptide fragment ions
including b, y, c, z ions with modifications, neutral losses, and
glycan-specific ions for EThcD spectrum annotation.

Author: Longping Fu
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set, Union
import json
import re

# Import glycan library for extended glycan support
from .glycan_library import (
    MONOSACCHARIDE_MASSES,
    OXONIUM_IONS_EXTENDED,
    OXONIUM_IONS_O_GLCNAC,
    OXONIUM_IONS_N_GLYCAN,
    OXONIUM_IONS_O_GLYCAN,
    GlycanComposition,
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    generate_y_ion_series,
    generate_n_glycan_y_ions,
    identify_glycan_from_mass,
)

# =============================================================================
# CONSTANTS: Monoisotopic Masses
# =============================================================================

# Amino acid residue masses (monoisotopic)
AA_MASSES = {
    'A': 71.03711,   # Alanine
    'R': 156.10111,  # Arginine
    'N': 114.04293,  # Asparagine
    'D': 115.02694,  # Aspartic acid
    'C': 103.00919,  # Cysteine
    'E': 129.04259,  # Glutamic acid
    'Q': 128.05858,  # Glutamine
    'G': 57.02146,   # Glycine
    'H': 137.05891,  # Histidine
    'I': 113.08406,  # Isoleucine
    'L': 113.08406,  # Leucine
    'K': 128.09496,  # Lysine
    'M': 131.04049,  # Methionine
    'F': 147.06841,  # Phenylalanine
    'P': 97.05276,   # Proline
    'S': 87.03203,   # Serine
    'T': 101.04768,  # Threonine
    'W': 186.07931,  # Tryptophan
    'Y': 163.06333,  # Tyrosine
    'V': 99.06841,   # Valine
}

# Modification masses
MOD_MASSES = {
    'TMT6plex': 229.1629,           # TMT6plex tag
    'HexNAc': 203.0794,             # N-acetylhexosamine (O-GlcNAc)
    'HexNAc_TMT': 528.2859,         # O-GlcNAc with TMT reporter
    'Carbamidomethyl': 57.02146,    # Cysteine alkylation
    'Oxidation': 15.9949,           # Methionine oxidation
}

# Common masses
PROTON = 1.007276
H2O = 18.010565
NH3 = 17.026549
CO = 27.994915
ELECTRON = 0.000549

# Neutral loss masses
NEUTRAL_LOSSES = {
    'H2O': 18.010565,
    'NH3': 17.026549,
    'CO2': 43.989829,         # Loss of CO2 (common for acidic residues)
    'HexNAc_TMT': 528.2859,   # Loss of glycan with TMT
    'HexNAc': 300.1308,       # Loss of glycan (oxonium form)
}

# Oxonium ions (glycan diagnostic ions) - B-type ions
# Basic set for O-GlcNAc (default for backwards compatibility)
OXONIUM_IONS = {
    'HexNAc_TMT': 529.2937,      # HexNAc with TMT oxonium
    'HexNAc+': 300.1308,         # HexNAc + H2O + H
    'HexNAc': 204.0864,          # HexNAc oxonium
    'HexNAc-H2O': 186.0760,      # HexNAc - H2O
    'HexNAc-2H2O': 168.0652,     # HexNAc - 2H2O
    'HexNAc_frag1': 144.0652,    # HexNAc fragment
    'HexNAc_frag2': 138.0546,    # HexNAc fragment
}

# Alias for extended oxonium ions (imported from glycan_library)
OXONIUM_IONS_ALL = OXONIUM_IONS_EXTENDED

# Glycan Y ion definitions (peptide + partial glycan)
# For O-GlcNAc with TMT: Y0 = peptide only, Y1 = peptide + HexNAc+TMT (intact)
# Mass values represent the LOSS from intact glycopeptide to get each Y ion
GLYCAN_Y_LOSSES = {
    'HexNAc_TMT': {
        # O-GlcNAc + TMT (528.2859 Da total)
        'Y0': 528.2859,           # Complete loss: peptide only
        'Y0-H2O': 528.2859 - 18.0106,  # Y0 with water loss from peptide
        'Y*': 203.0794,           # Loss of HexNAc only, TMT stays on peptide (rare)
    },
    'HexNAc': {
        # O-GlcNAc without TMT (203.0794 Da)
        'Y0': 203.0794,           # Complete loss: peptide only
        'Y0-H2O': 203.0794 - 18.0106,  # Y0 with water loss
    },
}

# Ion type mass adjustments (added to residue sum)
# For a fragment containing residues with total mass M:
ION_ADJUSTMENTS = {
    'b': PROTON,                           # +H+
    'y': H2O + PROTON,                     # +H2O +H+
    'c': NH3 + PROTON,                     # +NH3 +H+ (or +NH2 to peptide)
    'z': -NH3 + H2O + 2 * PROTON,          # z+H (even-electron, H rearrangement)
}

# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class Modification:
    """Represents a modification on a peptide."""
    position: int       # 0 = N-term, -1 = C-term, 1-based for residues
    residue: str        # 'N-term', 'C-term', or amino acid letter
    mass: float         # Modification mass
    name: str = ""      # Optional name


@dataclass
class TheoreticalIon:
    """Represents a theoretical fragment ion."""
    ion_type: str       # 'b', 'y', 'c', 'z', 'Y', 'oxonium', 'precursor'
    ion_number: int     # Ion number (e.g., b3 = 3)
    charge: int         # Charge state
    mz: float           # Theoretical m/z
    sequence: str       # Fragment sequence
    neutral_loss: str = ""  # Neutral loss type if any
    annotation: str = ""    # Full annotation string
    has_modification: bool = False  # Whether ion contains glycan/modification (marked with * in annotation)


@dataclass
class MatchedIon(TheoreticalIon):
    """Represents a matched ion with experimental data."""
    exp_mz: float = 0.0
    exp_intensity: float = 0.0
    mass_error_ppm: float = 0.0
    # Isotope-consistency flag (Phase 0 regression hardening). This is a
    # diagnostic FLAG, never a filter — the match stays in matched_ions.
    # Values:
    #   None: not evaluated
    #   'M1_absent': no M+1 isotope partner within tolerance (suspicious)
    #   'suspicious_low_M1': M+1 present but < 30% of averagine-expected ratio
    #   'suspicious_high_M1': M+1 > 300% of averagine-expected ratio (impossible)
    #   'ok': M+1 partner ratio consistent with averagine
    isotope_flag: Optional[str] = None
    # Phase 0b: precursor-envelope overlap audit fields. Set by
    # filter_precursor_envelope_overlap() when an ion is moved to the
    # excluded list. None on kept ions. See reviewer-approved Phase 0b
    # patch (April 2026): original-charge precursor isotope peaks frequently
    # coincide with fragment m/z (e.g. case 7's c3+N2 1+ landing on 4+
    # precursor M+2), producing spurious labels that the charge-reduced
    # filter — which only covers z < precursor_charge — cannot catch.
    exclusion_reason: Optional[str] = None
    exclusion_ppm: Optional[float] = None


# =============================================================================
# FRAGMENT CALCULATOR CLASS
# =============================================================================

class FragmentCalculator:
    """
    Calculates theoretical fragment ions for glycopeptides.

    Supports:
    - b, y ions (HCD)
    - c, z ions (ETD)
    - Y ions (precursor with glycan) - simple and extended series
    - Charge-reduced precursor species
    - Oxonium ions - basic and extended sets
    - Neutral losses (H2O, NH3, glycan)
    - Complex glycans (N-glycans, O-glycans)
    """

    def __init__(self,
                 peptide: str,
                 modifications: List[Dict],
                 precursor_charge: int,
                 max_fragment_charge: int = 2,
                 glycan_type: str = 'auto',
                 use_extended_oxonium: bool = False):
        """
        Initialize the fragment calculator.

        Args:
            peptide: Peptide sequence (single letter amino acids)
            modifications: List of dicts with 'position', 'residue', 'mass'
            precursor_charge: Precursor ion charge state
            max_fragment_charge: Maximum fragment ion charge (default 2)
            glycan_type: Type of glycan - 'auto', 'O-GlcNAc', 'O-GalNAc', 'N-glycan',
                        or a composition string like 'HexNAc2Hex5'
            use_extended_oxonium: Use extended oxonium ion library (for N-glycans)
        """
        self.peptide = peptide.upper()
        self.length = len(peptide)
        self.precursor_charge = precursor_charge
        self.max_fragment_charge = min(max_fragment_charge, precursor_charge)
        self.glycan_type_hint = glycan_type
        self.use_extended_oxonium = use_extended_oxonium

        # Parse modifications
        self.modifications = []
        self.mod_by_position = {}  # position -> Modification

        for mod in modifications:
            m = Modification(
                position=mod['position'],
                residue=mod.get('residue', ''),
                mass=mod['mass'],
                name=mod.get('name', '')
            )
            self.modifications.append(m)
            self.mod_by_position[m.position] = m

        # Calculate residue masses including modifications
        self._calculate_residue_masses()

        # Calculate precursor mass
        self.precursor_mass = self._calculate_precursor_mass()
        self.precursor_mz = (self.precursor_mass + precursor_charge * PROTON) / precursor_charge

    def _calculate_residue_masses(self):
        """Calculate mass of each residue including any modifications."""
        self.residue_masses = []

        for i, aa in enumerate(self.peptide):
            pos = i + 1  # 1-based position
            mass = AA_MASSES.get(aa, 0)

            # Add modification mass if present
            if pos in self.mod_by_position:
                mass += self.mod_by_position[pos].mass

            self.residue_masses.append(mass)

        # Store N-term and C-term modifications separately
        self.nterm_mod_mass = self.mod_by_position.get(0, Modification(0, '', 0)).mass
        self.cterm_mod_mass = self.mod_by_position.get(-1, Modification(-1, '', 0)).mass

    def _calculate_precursor_mass(self) -> float:
        """Calculate the neutral monoisotopic mass of the precursor."""
        # Sum of residues + H2O (terminal groups) + modifications
        mass = sum(self.residue_masses) + H2O
        mass += self.nterm_mod_mass
        mass += self.cterm_mod_mass
        return mass

    def _get_glycan_positions(self) -> List[int]:
        """Find all positions with glycan modifications.

        Recognizes glycans by:
        1. Name containing glycan monosaccharide names (HexNAc, Hex, NeuAc, Fuc, etc.)
        2. Known glycan masses (528.2859 HexNAc+TMT, 203.0794 HexNAc, 299.123 OPair)
        3. Modifications on S/T with mass > 150 Da (likely glycan)
        """
        glycan_name_parts = {'hex', 'hexnac', 'neuac', 'neugc', 'fuc',
                             'glyc', 'glcnac', 'galnac', 'sialyl'}
        known_glycan_masses = [528.2859, 203.0794, 299.123]
        positions = []

        for pos, mod in self.mod_by_position.items():
            if pos <= 0:  # Skip N-term/C-term
                continue

            # Check by name
            name_lower = mod.name.lower() if mod.name else ''
            if any(gn in name_lower for gn in glycan_name_parts):
                positions.append(pos)
                continue

            # Check by known mass
            if any(abs(mod.mass - m) < 0.5 for m in known_glycan_masses):
                positions.append(pos)
                continue

            # Check by mass on S/T (likely glycan if > 150 Da)
            if mod.residue in ('S', 'T') and mod.mass > 150:
                positions.append(pos)
                continue

        return sorted(positions)

    def _get_glycan_position(self) -> Optional[int]:
        """Find the first glycan modification position (backward compatible)."""
        positions = self._get_glycan_positions()
        return positions[0] if positions else None

    def _get_glycan_label(self) -> str:
        """Get the label string for the glycan modification type."""
        for pos, mod in self.mod_by_position.items():
            if abs(mod.mass - MOD_MASSES['HexNAc_TMT']) < 0.1:
                return '+HexNAc'
            elif abs(mod.mass - MOD_MASSES['HexNAc']) < 0.1:
                return '+HexNAc'
            elif abs(mod.mass - 299.123) < 0.5:
                return '+HexNAc'
        return '+Glyc'

    def _get_glycan_mass_for_ion(self, ion_type: str, ion_number: int) -> float:
        """Get total glycan mass covered by a fragment ion.

        For b/c ions, sums glycan masses at positions <= ion_number.
        For y/z ions, sums glycan masses at positions >= (length - ion_number + 1).
        """
        glycan_positions = self._get_glycan_positions()
        total = 0.0
        for pos in glycan_positions:
            if ion_type in ('b', 'c') and pos <= ion_number:
                total += self.mod_by_position[pos].mass
            elif ion_type in ('y', 'z') and pos >= (self.length - ion_number + 1):
                total += self.mod_by_position[pos].mass
        return total

    def calculate_b_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate b ions (N-terminal, HCD).
        b_n = sum of first n residues + N-term mod + H+

        For glycopeptides, generates both glycan-retaining and bare (glycan-lost)
        versions. Glycan-retaining ions provide localization evidence; bare ions
        match the dominant HCD fragmentation pathway (glycan lost before backbone cleavage).
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.nterm_mod_mass
        glycan_positions = self._get_glycan_positions()
        glycan_label = self._get_glycan_label()

        for i in range(self.length - 1):  # b1 to b(n-1)
            cumulative_mass += self.residue_masses[i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['b'] - PROTON  # Neutral mass

            ion_number = i + 1
            has_glycan = any(pos <= ion_number for pos in glycan_positions) if glycan_positions else False
            glycan_marker = glycan_label if has_glycan else ''

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='b',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[:i+1],
                    annotation=f"b{ion_number}{glycan_marker}{'⁺' * charge}",
                    has_modification=has_glycan
                )
                ions.append(ion)

        return ions

    def calculate_y_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate y ions (C-terminal, HCD).
        y_n = sum of last n residues + C-term mod + H2O + H+

        For glycopeptides, generates both glycan-retaining and bare versions.
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.cterm_mod_mass
        glycan_positions = self._get_glycan_positions()
        glycan_label = self._get_glycan_label()

        for i in range(self.length - 1):  # y1 to y(n-1)
            cumulative_mass += self.residue_masses[self.length - 1 - i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['y'] - PROTON  # Neutral mass

            # Check if this ion contains any glycan site
            # y_n covers positions (length-n+1) to length
            # Contains glycan if any glycan_pos >= (length - n + 1)
            ion_number = i + 1
            has_glycan = any(ion_number >= (self.length - pos + 1) for pos in glycan_positions) if glycan_positions else False
            glycan_marker = glycan_label if has_glycan else ''

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='y',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[-(i+1):],
                    annotation=f"y{ion_number}{glycan_marker}{'⁺' * charge}",
                    has_modification=has_glycan
                )
                ions.append(ion)

        return ions

    def calculate_c_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate c ions (N-terminal, ETD).
        c_n = b_n + NH3 = sum of first n residues + N-term mod + NH3 + H+

        For glycopeptides, generates both glycan-retaining and bare versions.
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.nterm_mod_mass
        glycan_positions = self._get_glycan_positions()
        glycan_label = self._get_glycan_label()

        for i in range(self.length - 1):  # c1 to c(n-1)
            cumulative_mass += self.residue_masses[i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['c'] - PROTON  # Neutral mass

            ion_number = i + 1
            has_glycan = any(pos <= ion_number for pos in glycan_positions) if glycan_positions else False
            glycan_marker = glycan_label if has_glycan else ''

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='c',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[:i+1],
                    annotation=f"c{ion_number}{glycan_marker}{'⁺' * charge}",
                    has_modification=has_glycan
                )
                ions.append(ion)

        return ions

    def calculate_z_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """
        Calculate z ions (C-terminal, ETD).
        z_n = y_n - NH3 + H (z+H even-electron ion, H rearrangement)

        For glycopeptides, generates both glycan-retaining and bare versions.
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        ions = []
        cumulative_mass = self.cterm_mod_mass
        glycan_positions = self._get_glycan_positions()
        glycan_label = self._get_glycan_label()

        for i in range(self.length - 1):  # z1 to z(n-1)
            cumulative_mass += self.residue_masses[self.length - 1 - i]
            ion_mass = cumulative_mass + ION_ADJUSTMENTS['z'] - PROTON  # Neutral mass

            ion_number = i + 1
            has_glycan = any(ion_number >= (self.length - pos + 1) for pos in glycan_positions) if glycan_positions else False
            glycan_marker = glycan_label if has_glycan else ''

            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='z',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[-(i+1):],
                    annotation=f"z{ion_number}{glycan_marker}{'⁺' * charge}",
                    has_modification=has_glycan
                )
                ions.append(ion)

        return ions

    def calculate_Y_ions(self, charges: List[int] = None,
                         extended_series: bool = False) -> List[TheoreticalIon]:
        """
        Calculate glycan Y ions (peptide + partial/no glycan).

        Y ions in glycopeptide MS represent the peptide backbone with
        varying amounts of glycan attached:
        - Y0: Peptide only (complete glycan loss)
        - Y1: Peptide + 1 monosaccharide (for complex glycans)
        - Y(intact): Full glycopeptide

        For O-GlcNAc (single monosaccharide):
        - Y0 = peptide backbone (glycan completely lost)
        - Y1 = intact glycopeptide

        For N-glycans (extended_series=True):
        - Y0, Y1, Y2, Y3, Y4, Y(core), Y(intact) + fucose variants

        Args:
            charges: List of charge states to generate
            extended_series: Generate full Y ion ladder for complex glycans
        """
        if charges is None:
            charges = list(range(1, self.precursor_charge + 1))

        ions = []

        # Determine glycan type and properties
        glycan_type = self._get_glycan_type()
        glycan_position = self._get_glycan_position()
        glycan_mass = self._get_glycan_mass()

        # Calculate peptide mass without glycan
        peptide_mass_no_glycan = self.precursor_mass - glycan_mass if glycan_mass else self.precursor_mass

        # Check if we should use extended Y ion series
        use_extended = extended_series or self.glycan_type_hint == 'N-glycan'

        if use_extended and glycan_mass:
            # Use extended Y ion series. The previous >400 gate excluded
            # Tn (203 Da), Core 1 (365 Da), sialyl-Tn (406 Da until decomposition
            # is precise), and other small biologically important O-glycans.
            ions.extend(self._calculate_extended_Y_ions(
                peptide_mass_no_glycan, glycan_mass, charges
            ))
        elif glycan_type and glycan_type in GLYCAN_Y_LOSSES:
            # Use simple Y ion series for O-GlcNAc
            y_losses = GLYCAN_Y_LOSSES[glycan_type]

            for y_name, loss_mass in y_losses.items():
                y_mass = self.precursor_mass - loss_mass

                for charge in charges:
                    mz = (y_mass + charge * PROTON) / charge

                    if 'Y0' in y_name:
                        ion_number = 0
                        annotation = f"{y_name} {charge}+"
                    else:
                        ion_number = 1
                        annotation = f"Y* {charge}+"

                    ion = TheoreticalIon(
                        ion_type='Y',
                        ion_number=ion_number,
                        charge=charge,
                        mz=mz,
                        sequence=self.peptide,
                        annotation=annotation
                    )
                    ions.append(ion)

        # Always add intact glycopeptide at different charge states. Only the
        # monoisotopic m/z is generated — the natural isotope envelope of the
        # surviving precursor in MS2 is not independent evidence and the
        # MSFragger/OPair searches operate at zero isotope error, so adding
        # +1iso/+2iso labels would mislead the reader into treating envelope
        # peaks as separate Y(intact) hits.
        for charge in charges:
            mz = (self.precursor_mass + charge * PROTON) / charge
            ion = TheoreticalIon(
                ion_type='Y',
                ion_number=1 if glycan_type else 0,
                charge=charge,
                mz=mz,
                sequence=self.peptide,
                annotation=f"Y(intact) {charge}+" if glycan_type else f"[M+{charge}H]{charge}+"
            )
            ions.append(ion)

        return ions

    @staticmethod
    def _decompose_glycan_mass(glycan_mass: float, tol: float = 0.05,
                               max_n: int = 12, max_h: int = 12,
                               max_a: int = 8, max_f: int = 6
                               ) -> Optional[Tuple[int, int, int, int]]:
        """Decompose a glycan total mass into (HexNAc, Hex, NeuAc, Fuc) counts.

        Returns the smallest-total-monosaccharide decomposition within
        ``tol`` Da, or None if no combination matches. The bounds are
        generous defaults for typical O- and N-glycans.
        """
        N = MONOSACCHARIDE_MASSES['HexNAc']
        H = MONOSACCHARIDE_MASSES['Hex']
        A = MONOSACCHARIDE_MASSES['NeuAc']
        F = MONOSACCHARIDE_MASSES['Fuc']
        best = None
        best_total = None
        for nN in range(max_n + 1):
            mN = nN * N
            if mN > glycan_mass + 1:
                break
            for nH in range(max_h + 1):
                mNH = mN + nH * H
                if mNH > glycan_mass + 1:
                    break
                for nA in range(max_a + 1):
                    mNHA = mNH + nA * A
                    if mNHA > glycan_mass + 1:
                        break
                    for nF in range(max_f + 1):
                        m = mNHA + nF * F
                        if m > glycan_mass + 1:
                            break
                        if abs(m - glycan_mass) < tol:
                            total = nN + nH + nA + nF
                            if best_total is None or total < best_total:
                                best = (nN, nH, nA, nF)
                                best_total = total
        return best

    def _calculate_extended_Y_ions(self, peptide_mass: float,
                                    glycan_mass: float,
                                    charges: List[int]) -> List[TheoreticalIon]:
        """Calculate extended Y ion series for any glycan composition.

        Strategy: decompose ``glycan_mass`` into (HexNAc, Hex, NeuAc, Fuc)
        counts, then enumerate every biologically-plausible sub-composition
        as a Y-ion (peptide + glycan stub). Sub-compositions are filtered
        to require at least one HexNAc as a reducing-end anchor (true for
        both O- and N-glycans).

        Falls back to a simple chitobiose ladder for masses > 800 Da that
        cannot be decomposed (e.g. with NeuGc, sulfate, phosphate).
        """
        ions = []
        N = MONOSACCHARIDE_MASSES['HexNAc']
        H = MONOSACCHARIDE_MASSES['Hex']
        A = MONOSACCHARIDE_MASSES['NeuAc']
        F = MONOSACCHARIDE_MASSES['Fuc']

        composition = self._decompose_glycan_mass(glycan_mass)
        y_ladder: Dict[str, float] = {'Y0': peptide_mass}

        if composition is not None:
            n_N, n_H, n_A, n_F = composition
            for i in range(n_N + 1):
                for j in range(n_H + 1):
                    for k in range(n_A + 1):
                        for l in range(n_F + 1):
                            if (i, j, k, l) == (0, 0, 0, 0):
                                continue  # That's Y0, already added
                            if (i, j, k, l) == (n_N, n_H, n_A, n_F):
                                continue  # Y(intact), added by calculate_Y_ions
                            # Require at least one HexNAc anchor
                            # (reducing-end of O- and N-glycans is HexNAc)
                            if i == 0 and (j > 0 or k > 0 or l > 0):
                                continue
                            stub_mass = i * N + j * H + k * A + l * F
                            parts = []
                            if i: parts.append(f'N{i}' if i > 1 else 'N')
                            if j: parts.append(f'H{j}' if j > 1 else 'H')
                            if k: parts.append(f'A{k}' if k > 1 else 'A')
                            if l: parts.append(f'F{l}' if l > 1 else 'F')
                            label = 'Y0+' + ''.join(parts)
                            y_ladder[label] = peptide_mass + stub_mass

        elif glycan_mass > 800:
            # Mass not decomposable — fall back to chitobiose ladder
            y_ladder['Y0+N'] = peptide_mass + N
            y_ladder['Y0+N2'] = peptide_mass + 2 * N
            y_ladder['Y0+N2H'] = peptide_mass + 2 * N + H
            y_ladder['Y0+N2H2'] = peptide_mass + 2 * N + 2 * H
            y_ladder['Y0+N2H3'] = peptide_mass + 2 * N + 3 * H
            if glycan_mass > 1200:
                y_ladder['Y0+N2H3F'] = peptide_mass + 2 * N + 3 * H + F
        else:
            # Tiny / unknown glycan mass — at least add the HexNAc stub
            y_ladder['Y0+N'] = peptide_mass + N

        # Generate TheoreticalIon objects
        for y_name, y_mass in y_ladder.items():
            if y_mass > self.precursor_mass + 1:
                continue

            ion_number = 0 if y_name == 'Y0' else 1

            for charge in charges:
                mz = (y_mass + charge * PROTON) / charge
                ion = TheoreticalIon(
                    ion_type='Y',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide,
                    annotation=f"{y_name} {charge}+"
                )
                ions.append(ion)

        return ions

    def _get_glycan_mass(self) -> float:
        """Get the total mass of glycan modification(s).

        Delegates detection to ``_get_glycan_positions`` so that the
        positions, masses, and types are all derived from a single source
        of truth. This catches medium-mass O-glycans (e.g. Core 1 365 Da,
        sialyl-Tn 406 Da) that the previous implementation silently
        returned 0 for, which broke the extended Y-ladder generation.
        """
        glycan_positions = self._get_glycan_positions()
        return sum(self.mod_by_position[p].mass for p in glycan_positions)

    def _get_glycan_type(self) -> Optional[str]:
        """
        Determine the type of glycan modification.

        Returns:
            Glycan type string: 'HexNAc_TMT', 'HexNAc', 'O-glycan', 'N-glycan', or None

        Important: a modification on S/T is treated as O-glycan regardless of
        mass. The previous implementation classified any glycan > 800 Da as
        'N-glycan' even when the modification site was S/T, which then caused
        the wrong (high-mannose) oxonium ion set to be used for annotation.
        """
        # If user specified glycan type, use that
        if self.glycan_type_hint and self.glycan_type_hint != 'auto':
            return self.glycan_type_hint

        # Multi-site glycopeptides: classify by the LARGEST glycan stub, not
        # the first one in iteration order. A peptide with stubs [203, 1022,
        # 860] is a complex multi-site O-glycan, not "HexNAc/O-GlcNAc" — but
        # the previous code returned on the first 203 match and selected the
        # simplified O-GlcNAc oxonium set, which excludes NeuAc/Fuc diagnostics.
        glycan_mods = [m for m in self.mod_by_position.values()
                       if (m.residue in ('S', 'T') and m.mass > 150)
                       or (m.residue == 'N' and m.mass > 800)
                       or m.mass > 200]
        if not glycan_mods:
            return None

        largest = max(glycan_mods, key=lambda m: m.mass)
        # Single-site O-GlcNAc detection by exact mass on the largest stub
        if abs(largest.mass - MOD_MASSES['HexNAc_TMT']) < 0.1:
            return 'HexNAc_TMT'
        if abs(largest.mass - MOD_MASSES['HexNAc']) < 0.1:
            return 'HexNAc'
        if abs(largest.mass - 299.123) < 0.5:
            return 'HexNAc'
        # If any glycosite is on S/T with mass > 150, it's a mucin-type O-glycan
        if any(m.residue in ('S', 'T') and m.mass > 150 for m in glycan_mods):
            return 'O-glycan'
        # N-glycan if largest mod is on N residue
        if largest.residue == 'N' and largest.mass > 800:
            return 'N-glycan'
        # Last-resort fallbacks
        if largest.mass > 800:
            return 'N-glycan'
        if largest.mass > 350:
            return 'O-glycan'
        return None

    def calculate_charge_reduced_precursor(self) -> List[TheoreticalIon]:
        """Calculate charge-reduced precursor species from ETD.

        Includes:
        - Charge-reduced radical [M+nH](n-1)+• at each lower charge state
        - Neutral losses from charge-reduced: -H2O, -NH3, -H2 (2.016 Da)
        - Isotope peaks (M+1, M+2) for each species
        """
        ions = []
        neutral = self.precursor_mass

        losses = [
            (0.0, ''),
            (18.010565, '-H2O'),
            (17.026549, '-NH3'),
            (2.01565, '-H2'),
        ]

        precursor_protons = self.precursor_charge * PROTON

        for z in range(1, self.precursor_charge):
            for loss_mass, loss_label in losses:
                for iso in range(3):
                    mz = (neutral + precursor_protons - loss_mass) / z + iso * 1.003355 / z
                    iso_label = f"+{iso}" if iso > 0 else ""
                    if loss_label:
                        ann = f"[M{loss_label}+{self.precursor_charge}H]{z}+\u2022{iso_label}"
                    else:
                        ann = f"[M+{self.precursor_charge}H]{z}+\u2022{iso_label}"
                    ions.append(TheoreticalIon(
                        ion_type='precursor', ion_number=iso, charge=z,
                        mz=mz, sequence=self.peptide, annotation=ann))

        # Losses at original charge
        z = self.precursor_charge
        for loss_mass, loss_label in losses[1:]:
            for iso in range(3):
                mz = (neutral - loss_mass + z * PROTON) / z + iso * 1.003355 / z
                iso_label = f"+{iso}" if iso > 0 else ""
                ann = f"[M{loss_label}+{z}H]{z}+{iso_label}"
                ions.append(TheoreticalIon(
                    ion_type='precursor', ion_number=iso, charge=z,
                    mz=mz, sequence=self.peptide, annotation=ann))

        return ions

    def calculate_precursor_isotopes(self, n_isotopes: int = 4) -> List[TheoreticalIon]:
        """Calculate precursor isotope peaks."""
        ions = []
        isotope_spacing = 1.003355 / self.precursor_charge

        for i in range(n_isotopes):
            mz = self.precursor_mz + i * isotope_spacing
            ion = TheoreticalIon(
                ion_type='precursor',
                ion_number=i,
                charge=self.precursor_charge,
                mz=mz,
                sequence=self.peptide,
                annotation=f"[M+{self.precursor_charge}H]{self.precursor_charge}+ iso{i}"
            )
            ions.append(ion)

        return ions

    def calculate_oxonium_ions(self, use_extended: bool = None) -> List[TheoreticalIon]:
        """
        Calculate glycan oxonium (diagnostic) ions.

        Args:
            use_extended: Use extended oxonium ion library. If None, uses
                         self.use_extended_oxonium setting.

        Returns:
            List of oxonium ion TheoreticalIon objects
        """
        ions = []

        if use_extended is None:
            use_extended = self.use_extended_oxonium

        # Select oxonium ion set based on glycan type and settings
        if use_extended:
            oxonium_set = OXONIUM_IONS_EXTENDED
        else:
            # Auto-select based on detected glycan type. Note: any glycan
            # whose modification site is on S/T should be treated as O-glycan
            # regardless of mass; high-mannose-style ions
            # (HexNAc-2Hex/HexNAc-3Hex/HexNAc-4Hex) do not occur in
            # mucin-type O-glycan fragmentation.
            glycan_type = self._get_glycan_type()
            if glycan_type in ['HexNAc_TMT', 'HexNAc', 'O-GlcNAc', 'O-GalNAc']:
                oxonium_set = OXONIUM_IONS_O_GLCNAC
            elif glycan_type == 'O-glycan':
                oxonium_set = OXONIUM_IONS_O_GLYCAN
            elif glycan_type == 'N-glycan':
                oxonium_set = OXONIUM_IONS_N_GLYCAN
            else:
                oxonium_set = OXONIUM_IONS

        for name, mz in oxonium_set.items():
            ion = TheoreticalIon(
                ion_type='oxonium',
                ion_number=0,
                charge=1,
                mz=mz,
                sequence='',
                annotation=name
            )
            ions.append(ion)

        return ions

    # Residues that can produce specific neutral losses
    _H2O_RESIDUES = set('STDE')   # Ser, Thr (hydroxyl), Asp, Glu (carboxyl)
    _NH3_RESIDUES = set('RKNQ')   # Arg, Lys (amine), Asn, Gln (amide)

    def _get_fragment_sequence(self, ion_type: str, ion_number: int) -> str:
        """Get the amino acid sequence covered by a fragment ion."""
        if ion_type in ('b', 'c'):
            return self.peptide[:ion_number]
        elif ion_type in ('y', 'z'):
            return self.peptide[self.length - ion_number:]
        return ''

    def calculate_neutral_loss_ions(self,
                                    base_ions: List[TheoreticalIon],
                                    loss_types: List[str] = None) -> List[TheoreticalIon]:
        """
        Calculate neutral loss ions from base fragment ions.
        Only applies losses when the fragment contains applicable residues:
          H2O loss: requires S, T, D, or E in the fragment
          NH3 loss: requires R, K, N, or Q in the fragment

        Args:
            base_ions: List of base fragment ions
            loss_types: List of neutral loss types ('H2O', 'NH3', 'HexNAc_TMT', 'HexNAc')
        """
        if loss_types is None:
            loss_types = ['H2O', 'NH3']

        glycan_pos = self._get_glycan_position()
        ions = []

        for base_ion in base_ions:
            # Get the sequence of this fragment for residue-specific checks
            frag_seq = self._get_fragment_sequence(base_ion.ion_type, base_ion.ion_number)

            for loss_type in loss_types:
                loss_mass = NEUTRAL_LOSSES.get(loss_type, 0)
                if loss_mass == 0:
                    continue

                # Residue-specific filtering for H2O and NH3 losses
                if loss_type == 'H2O':
                    if not any(aa in self._H2O_RESIDUES for aa in frag_seq):
                        continue
                elif loss_type == 'NH3':
                    # Skip NH3 loss for c ions: c_n - NH3 = b_n (identical mass)
                    if base_ion.ion_type == 'c':
                        continue
                    if not any(aa in self._NH3_RESIDUES for aa in frag_seq):
                        continue

                # For glycan loss, only apply if the fragment contains the glycan
                if loss_type in ['HexNAc_TMT', 'HexNAc'] and glycan_pos:
                    if base_ion.ion_type in ['b', 'c']:
                        if base_ion.ion_number < glycan_pos:
                            continue
                    elif base_ion.ion_type in ['y', 'z']:
                        if base_ion.ion_number < (self.length - glycan_pos + 1):
                            continue

                # Calculate neutral loss m/z
                new_mz = base_ion.mz - loss_mass / base_ion.charge

                if new_mz > 0:
                    ion = TheoreticalIon(
                        ion_type=base_ion.ion_type,
                        ion_number=base_ion.ion_number,
                        charge=base_ion.charge,
                        mz=new_mz,
                        sequence=base_ion.sequence,
                        neutral_loss=loss_type,
                        annotation=f"{base_ion.annotation}-{loss_type}",
                        has_modification=base_ion.has_modification
                    )
                    ions.append(ion)

        return ions

    def calculate_all_ions(self,
                          include_neutral_losses: bool = True,
                          neutral_loss_types: List[str] = None,
                          extended_y_series: bool = False) -> Dict[str, List[TheoreticalIon]]:
        """
        Calculate all theoretical ions for the peptide (including c/z for EThcD).

        Note: Glycan neutral losses are NOT included by default because they
        don't provide localization information.

        Returns dict organized by ion type.
        """
        if neutral_loss_types is None:
            # Only H2O and NH3 losses - NO glycan losses
            neutral_loss_types = ['H2O', 'NH3']

        result = {
            'b': self.calculate_b_ions(),
            'y': self.calculate_y_ions(),
            'c': self.calculate_c_ions(),
            'z': self.calculate_z_ions(),
            'Y': self.calculate_Y_ions(extended_series=extended_y_series),
            'precursor': self.calculate_precursor_isotopes() + self.calculate_charge_reduced_precursor(),
            'oxonium': self.calculate_oxonium_ions(),
        }

        if include_neutral_losses:
            # Add neutral losses for backbone ions
            for ion_type in ['b', 'y', 'c', 'z']:
                nl_ions = self.calculate_neutral_loss_ions(result[ion_type], neutral_loss_types)
                result[f'{ion_type}_NL'] = nl_ions

        return result

    def calculate_glycan_loss_by_ions(self, charges: List[int] = None) -> List[TheoreticalIon]:
        """Bare (glycan-lost) b/y ions for fragments crossing glycosites.

        Under HCD, glycopeptides often undergo complete glycan loss *before*
        backbone cleavage, producing b/y ions whose m/z matches the naked
        peptide value even when the fragment would nominally span a glycosite.
        These ions are key peptide-mass evidence (same principle as Y0) and
        should appear in the theoretical-ion set so the annotator can label
        them — otherwise only glycan-bearing b/y at out-of-range m/z are
        computed, and the glycan-loss peaks stay unlabeled.

        Only ions that cross at least one glycosite are returned (fragments
        that don't cross a glycosite are identical to normal calculate_b_ions
        / calculate_y_ions output, so emitting them again would duplicate).
        """
        if charges is None:
            charges = list(range(1, self.max_fragment_charge + 1))

        glycan_positions = self._get_glycan_positions()
        if not glycan_positions:
            return []

        # Build a "naked" residue-mass array: for glycosylated S/T residues,
        # subtract the glycan mass.
        naked_residues = list(self.residue_masses)
        for pos in glycan_positions:
            mod = self.mod_by_position.get(pos)
            if mod is not None:
                naked_residues[pos - 1] -= mod.mass

        ions: List[TheoreticalIon] = []

        # Bare b-ions: generate for b_n where n >= min(glycan_positions)
        # (ions below the first glycosite are already bare in calculate_b_ions).
        cum = self.nterm_mod_mass
        for i in range(self.length - 1):
            cum += naked_residues[i]
            ion_number = i + 1
            if ion_number < min(glycan_positions):
                continue
            ion_mass = cum + ION_ADJUSTMENTS['b'] - PROTON
            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ions.append(TheoreticalIon(
                    ion_type='b',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[:i+1],
                    annotation=f"b{ion_number}°{'⁺' * charge}",
                    has_modification=False,
                    neutral_loss='glyc_full',
                ))

        # Bare y-ions: mirror calculate_y_ions exactly, but with naked residues.
        # y_n covers the last n residues; contains glycan if any glycan_pos
        # is >= (length - n + 1). Skip n where no glycosite is crossed
        # (those match calculate_y_ions output already).
        cum = self.cterm_mod_mass
        for i in range(self.length - 1):
            cum += naked_residues[self.length - 1 - i]
            ion_number = i + 1
            has_glycan = any(ion_number >= (self.length - pos + 1)
                              for pos in glycan_positions)
            if not has_glycan:
                continue
            ion_mass = cum + ION_ADJUSTMENTS['y'] - PROTON
            for charge in charges:
                mz = (ion_mass + charge * PROTON) / charge
                ions.append(TheoreticalIon(
                    ion_type='y',
                    ion_number=ion_number,
                    charge=charge,
                    mz=mz,
                    sequence=self.peptide[-(i+1):],
                    annotation=f"y{ion_number}°{'⁺' * charge}",
                    has_modification=False,
                    neutral_loss='glyc_full',
                ))

        return ions

    def calculate_hcd_ions(self,
                          include_neutral_losses: bool = True,
                          neutral_loss_types: List[str] = None,
                          extended_y_series: bool = False,
                          include_glycan_loss_by: bool = True) -> Dict[str, List[TheoreticalIon]]:
        """
        Calculate theoretical ions for HCD spectra only.

        Includes: b, y, Y (intact glycopeptide), precursor, oxonium
        Excludes: c, z (ETD-specific ions)

        Note: For glycopeptides, glycan neutral losses (HexNAc_TMT, HexNAc) are NOT
        included because they result in bare peptide ions that don't provide
        localization information. Only H2O and NH3 losses are included.

        include_glycan_loss_by: when True (default), also emit bare b/y ions
        for fragments crossing glycosites — these are the "glycan-loss" HCD
        products that would otherwise be unmatched because the standard
        calculate_b_ions / calculate_y_ions compute theoretical m/z with the
        full glycan attached. See `calculate_glycan_loss_by_ions`.

        Returns dict organized by ion type.
        """
        if neutral_loss_types is None:
            # Only H2O and NH3 losses - NO glycan losses for localization analysis
            neutral_loss_types = ['H2O', 'NH3']

        result = {
            'b': self.calculate_b_ions(),
            'y': self.calculate_y_ions(),
            'Y': self.calculate_Y_ions(extended_series=extended_y_series),
            'precursor': self.calculate_precursor_isotopes() + self.calculate_charge_reduced_precursor(),
            'oxonium': self.calculate_oxonium_ions(),
        }

        if include_glycan_loss_by:
            result['by_glycloss'] = self.calculate_glycan_loss_by_ions()

        if include_neutral_losses:
            # Add neutral losses for backbone ions (b and y only, no c/z)
            for ion_type in ['b', 'y']:
                nl_ions = self.calculate_neutral_loss_ions(result[ion_type], neutral_loss_types)
                result[f'{ion_type}_NL'] = nl_ions

        return result

    def get_hcd_ions_flat(self, **kwargs) -> List[TheoreticalIon]:
        """Get HCD-appropriate ions as a flat list (no c/z ions)."""
        all_ions = self.calculate_hcd_ions(**kwargs)
        flat_list = []
        for ion_list in all_ions.values():
            flat_list.extend(ion_list)
        return flat_list

    def get_all_ions_flat(self, **kwargs) -> List[TheoreticalIon]:
        """Get all ions as a flat list."""
        all_ions = self.calculate_all_ions(**kwargs)
        flat_list = []
        for ion_list in all_ions.values():
            flat_list.extend(ion_list)
        return flat_list


# =============================================================================
# DEISOTOPING
# =============================================================================

ISOTOPE_SPACING = 1.003355


def deisotope(exp_mz: np.ndarray,
              exp_intensity: np.ndarray,
              max_charge: int = 4,
              tolerance_ppm: float = 15.0,
              min_isotope_ratio: float = 0.01,
              max_isotope_ratio: float = 2.0,
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove isotope peaks from a centroided spectrum, keeping monoisotopic peaks.

    .. deprecated::
        Use ``mzml_utils.deisotope()`` instead, which uses averagine-scored
        cluster detection with subcluster validation.
    """
    import warnings
    warnings.warn(
        "Use mzml_utils.deisotope() instead (averagine-scored deisotoping)",
        DeprecationWarning, stacklevel=2,
    )
    if len(exp_mz) == 0:
        return exp_mz, exp_intensity

    is_isotope = np.zeros(len(exp_mz), dtype=bool)

    for i in range(len(exp_mz)):
        if is_isotope[i]:
            continue
        for z in range(1, max_charge + 1):
            spacing = ISOTOPE_SPACING / z
            target_mz = exp_mz[i] + spacing
            tol = target_mz * tolerance_ppm / 1e6
            for j in range(i + 1, len(exp_mz)):
                if exp_mz[j] > target_mz + tol:
                    break
                if abs(exp_mz[j] - target_mz) <= tol and not is_isotope[j]:
                    ratio = exp_intensity[j] / exp_intensity[i] if exp_intensity[i] > 0 else 0
                    if min_isotope_ratio <= ratio <= max_isotope_ratio:
                        is_isotope[j] = True
                        target_m2 = exp_mz[i] + 2 * spacing
                        tol_m2 = target_m2 * tolerance_ppm / 1e6
                        for k in range(j + 1, len(exp_mz)):
                            if exp_mz[k] > target_m2 + tol_m2:
                                break
                            if abs(exp_mz[k] - target_m2) <= tol_m2 and not is_isotope[k]:
                                if exp_intensity[k] / exp_intensity[i] <= max_isotope_ratio:
                                    is_isotope[k] = True
                                    break
                    break

    keep = ~is_isotope
    return exp_mz[keep], exp_intensity[keep]


# =============================================================================
# CHARGE-REDUCED PRECURSOR EXCLUSION (EThcD)
# =============================================================================

def build_charge_reduced_exclusions(
    precursor_neutral_mass: float,
    charge: int,
    glycan_mass: float,
    n_isotopes: int = 4,
) -> List[float]:
    """Build exclusion list for charge-reduced precursor species (ETnoD)."""
    exclusions = []
    neutral_losses = [0.0, 203.07937, 162.05282, 291.09542]
    for reduced_z in range(1, charge):
        for nl_mass in neutral_losses:
            if nl_mass > glycan_mass + 0.01:
                continue
            for iso in range(n_isotopes):
                cr_mz = (precursor_neutral_mass - nl_mass
                         + reduced_z * PROTON
                         + iso * ISOTOPE_SPACING) / reduced_z
                exclusions.append(cr_mz)
    return exclusions


def filter_charge_reduced(
    matched_ions: List['MatchedIon'],
    exclusion_mzs: List[float],
    exclusion_ppm: float = 15.0,
) -> List['MatchedIon']:
    """Remove matched c/z ions that overlap with charge-reduced precursor species."""
    if not exclusion_mzs:
        return matched_ions
    exclusion_arr = np.array(exclusion_mzs)
    filtered = []
    for m in matched_ions:
        if m.ion_type in ('c', 'z'):
            ppm_diffs = np.abs(m.exp_mz - exclusion_arr) / exclusion_arr * 1e6
            if np.min(ppm_diffs) < exclusion_ppm:
                continue
        filtered.append(m)
    return filtered


# =============================================================================
# ORIGINAL-CHARGE PRECURSOR ENVELOPE EXCLUSION (Phase 0b)
# =============================================================================
#
# Reviewer-approved follow-up to Phase 0 (April 2026). build_charge_reduced
# _exclusions() only covers z < precursor_charge, but spurious fragment labels
# also land on the ORIGINAL charge precursor isotope envelope. Example: case 7
# has a 4+ precursor at 986.92; M+2 lands at 987.42 (986.92 + 2 * 0.251). A
# theoretical c3+N2 1+ at 987.45 ppm-matches this peak and gets labelled even
# though the peak is really the surviving precursor's 3rd isotopologue.

def build_precursor_envelope_exclusions(
    precursor_mz: float,
    precursor_charge: int,
    n_isotopes_below: int = 2,
    n_isotopes_above: int = 8,
) -> List[float]:
    """Discrete-centroid exclusion list for the ORIGINAL-charge precursor isotope envelope.

    Complements build_charge_reduced_exclusions() which only covers z < charge.
    The original-charge isotope peaks can coincide with fragment m/z (e.g. c3+N2 1+
    landing on 4+ precursor M+2), producing spurious matches. Returns discrete
    centroid m/z values; callers apply a narrow ppm tolerance around each centroid.

    Defaults cover M-2 through M+8 to handle large glycopeptides (4-6 kDa) where
    the averagine envelope is wide and MIPS offsets can push observed below mono.
    """
    return [
        precursor_mz + k * ISOTOPE_SPACING / precursor_charge
        for k in range(-n_isotopes_below, n_isotopes_above + 1)
    ]


def filter_precursor_envelope_overlap(
    matched_ions: List['MatchedIon'],
    exclusion_mzs: List[float],
    precursor_mz: float,
    precursor_charge: int,
    exclusion_ppm: float = 30.0,
) -> Tuple[List['MatchedIon'], List['MatchedIon']]:
    """Remove matches whose m/z falls within exclusion_ppm of precursor-envelope
    centroids. Returns (filtered, excluded) — excluded matches carry
    ``exclusion_reason='original_precursor_envelope_overlap'`` and
    ``exclusion_ppm`` attributes for audit.

    Preserved:
      - oxonium ions (low-mass, unrelated to precursor envelope range)
      - Y(intact) at precursor_charge (legitimate surviving precursor);
        identified by annotation starting with 'Y(intact)' AND m/z within
        exclusion_ppm of the observed precursor mono m/z.

    Note: all other Y-type ions (Y0, Y-intermediates) currently share
    ion_number=1 with Y(intact), so we MUST NOT distinguish by ion_number alone.
    Use annotation string + m/z proximity to the selected precursor m/z.
    """
    if not exclusion_mzs:
        return list(matched_ions), []
    exclusion_arr = np.array(exclusion_mzs)
    filtered = []
    excluded = []
    for m in matched_ions:
        if m.ion_type == "oxonium":
            filtered.append(m)
            continue
        # Recognize Y(intact) at precursor_charge by annotation + m/z.
        # Reviewer correction: do NOT distinguish by ion_number alone because
        # Y0 / Y-intermediates / Y(intact) all share ion_number=1 in this package.
        is_intact_y = (
            m.ion_type == "Y"
            and m.charge == precursor_charge
            and m.annotation.startswith("Y(intact)")
            and abs(m.exp_mz - precursor_mz) / precursor_mz * 1e6 < exclusion_ppm
        )
        if is_intact_y:
            filtered.append(m)
            continue
        ppm_diffs = np.abs(m.exp_mz - exclusion_arr) / exclusion_arr * 1e6
        min_ppm = float(np.min(ppm_diffs))
        if min_ppm < exclusion_ppm:
            m.exclusion_reason = "original_precursor_envelope_overlap"
            m.exclusion_ppm = min_ppm
            excluded.append(m)
            continue
        filtered.append(m)
    return filtered, excluded


# =============================================================================
# S/N CALCULATION
# =============================================================================

def load_noise_cache(noise_mzml_path: str) -> Dict[int, Dict[str, np.ndarray]]:
    """Load noise profiles from ThermoRawFileParser -N mzML file."""
    from pyteomics import mzml
    cache = {}
    reader = mzml.MzML(str(noise_mzml_path))
    for spec in reader:
        scan_id = spec.get('id', '')
        if 'scan=' in scan_id:
            scan_num = int(scan_id.split('scan=')[-1])
        else:
            continue
        if 'sampled noise m/z array' in spec:
            cache[scan_num] = {
                'noise_mz': spec['sampled noise m/z array'],
                'noise_int': spec['sampled noise intensity array'],
            }
    reader.close()
    return cache


def get_sn_array(exp_mz: np.ndarray,
                 exp_intensity: np.ndarray,
                 noise_cache: Optional[Dict] = None,
                 scan_num: Optional[int] = None,
                 ) -> np.ndarray:
    """Calculate S/N per peak. Uses instrument noise if available, else median."""
    if len(exp_mz) == 0:
        return np.array([])
    if noise_cache is not None and scan_num is not None:
        noise_data = noise_cache.get(scan_num)
        if noise_data is not None:
            noise_at_peaks = np.interp(exp_mz, noise_data['noise_mz'],
                                        noise_data['noise_int'])
            noise_at_peaks = np.maximum(noise_at_peaks, 1.0)
            return exp_intensity / noise_at_peaks
    noise_level = max(np.median(exp_intensity), 1.0)
    return exp_intensity / noise_level


# =============================================================================
# PEAK MATCHING
# =============================================================================

def match_peaks(theoretical_ions: List[TheoreticalIon],
                exp_mz: np.ndarray,
                exp_intensity: np.ndarray,
                tolerance_ppm: float = 20.0,
                match_isotopes: bool = True,
                max_isotope: int = 2) -> List[MatchedIon]:
    """
    Match experimental peaks to theoretical ions, including isotope peaks.

    Args:
        theoretical_ions: List of theoretical ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        tolerance_ppm: Mass tolerance in ppm
        match_isotopes: Whether to also match M+1, M+2 isotope peaks
        max_isotope: Maximum isotope offset to check (default 2 for M+0, M+1, M+2)

    Returns:
        List of matched ions with experimental data
    """
    matched = []

    # Isotope mass spacing (C13 - C12)
    ISOTOPE_SPACING = 1.003355

    # Sort theoretical ions by priority:
    # 1. Y ions (glycan) first - these are important diagnostic ions
    # 2. Base ions before neutral loss ions
    # 3. Then by m/z
    def ion_sort_key(ion):
        # Priority: Y ions get highest priority (0), then base ions (1), then NL ions (2)
        if ion.ion_type == 'Y':
            type_priority = 0
        elif not ion.neutral_loss:
            type_priority = 1
        else:
            type_priority = 2
        return (type_priority, ion.mz)

    sorted_ions = sorted(theoretical_ions, key=ion_sort_key)

    for ion in sorted_ions:
        theo_mz = ion.mz
        found_match = False

        # Try monoisotopic first, then isotopes if enabled
        isotope_offsets = [0]
        if match_isotopes:
            isotope_offsets = list(range(max_isotope + 1))  # [0, 1, 2]

        for iso_offset in isotope_offsets:
            if found_match:
                break

            # Calculate m/z with isotope offset (divided by charge)
            target_mz = theo_mz + (iso_offset * ISOTOPE_SPACING / ion.charge)
            tol = target_mz * tolerance_ppm / 1e6

            # Find peaks within tolerance
            matches = np.where(np.abs(exp_mz - target_mz) <= tol)[0]

            if len(matches) > 0:
                # Find the closest match (allow multiple ions to match same peak)
                best_idx = matches[np.argmin(np.abs(exp_mz[matches] - target_mz))]
                error_ppm = (exp_mz[best_idx] - target_mz) / target_mz * 1e6

                annotation = ion.annotation
                if iso_offset > 0:
                    annotation = f"{ion.annotation}+{iso_offset}"

                matched_ion = MatchedIon(
                    ion_type=ion.ion_type,
                    ion_number=ion.ion_number,
                    charge=ion.charge,
                    mz=ion.mz,
                    sequence=ion.sequence,
                    neutral_loss=ion.neutral_loss,
                    annotation=annotation,
                    has_modification=ion.has_modification,
                    exp_mz=exp_mz[best_idx],
                    exp_intensity=exp_intensity[best_idx],
                    mass_error_ppm=error_ppm
                )
                matched.append(matched_ion)
                found_match = True
                break

    return matched


# =============================================================================
# FALSE MATCH RATE CALCULATION
# =============================================================================

@dataclass
class FalseMatchRate:
    """Results of false match rate estimation."""
    fmr_peaks: float          # False match rate based on peak count
    fmr_intensity: float      # False match rate based on intensity
    matched_peaks: int        # Number of matched peaks in true spectrum
    matched_intensity: float  # Sum of matched intensity in true spectrum
    avg_random_peaks: float   # Average matched peaks in shifted spectra
    avg_random_intensity: float  # Average matched intensity in shifted spectra
    n_shifts: int             # Number of spectrum shifts used


def calculate_false_match_rate(
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    tolerance_ppm: float = 20.0,
    shift_range: float = 25.0,
    shift_step: float = 1.0
) -> FalseMatchRate:
    """
    Calculate false match rate using spectrum shifting method.

    This implements the method from Schulte et al. (Anal. Chem. 2025):
    The spectrum is shifted by π ± shift_range Th in shift_step increments.
    The π offset prevents false matches from isotope patterns.

    Args:
        theoretical_ions: List of theoretical fragment ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        tolerance_ppm: Mass tolerance for matching (default 20 ppm)
        shift_range: Range of shifts in Th (default ±25 Th)
        shift_step: Step size for shifts (default 1 Th)

    Returns:
        FalseMatchRate object with peak and intensity-based FMR
    """
    # First, match the true (unshifted) spectrum
    true_matched = match_peaks(
        theoretical_ions, exp_mz, exp_intensity,
        tolerance_ppm, match_isotopes=False  # Disable isotope matching for FMR
    )

    true_peak_count = len(true_matched)
    true_intensity_sum = sum(ion.exp_intensity for ion in true_matched)

    # If no matches in true spectrum, return zero FMR
    if true_peak_count == 0:
        return FalseMatchRate(
            fmr_peaks=0.0,
            fmr_intensity=0.0,
            matched_peaks=0,
            matched_intensity=0.0,
            avg_random_peaks=0.0,
            avg_random_intensity=0.0,
            n_shifts=0
        )

    # Generate shift offsets: π ± shift_range in shift_step increments
    # π is added to prevent false matches from isotope patterns (spacing ~1 Da)
    pi_offset = np.pi
    shifts = np.arange(-shift_range, shift_range + shift_step, shift_step) + pi_offset

    # Track matches for each shifted spectrum
    shifted_peak_counts = []
    shifted_intensity_sums = []

    for shift in shifts:
        # Shift the experimental m/z values
        shifted_mz = exp_mz + shift

        # Match against theoretical ions
        shifted_matched = match_peaks(
            theoretical_ions, shifted_mz, exp_intensity,
            tolerance_ppm, match_isotopes=False
        )

        shifted_peak_counts.append(len(shifted_matched))
        shifted_intensity_sums.append(sum(ion.exp_intensity for ion in shifted_matched))

    # Calculate average matches across all shifts
    avg_random_peaks = np.mean(shifted_peak_counts)
    avg_random_intensity = np.mean(shifted_intensity_sums)

    # Calculate false match rates
    fmr_peaks = avg_random_peaks / true_peak_count if true_peak_count > 0 else 0.0
    fmr_intensity = avg_random_intensity / true_intensity_sum if true_intensity_sum > 0 else 0.0

    return FalseMatchRate(
        fmr_peaks=fmr_peaks,
        fmr_intensity=fmr_intensity,
        matched_peaks=true_peak_count,
        matched_intensity=true_intensity_sum,
        avg_random_peaks=avg_random_peaks,
        avg_random_intensity=avg_random_intensity,
        n_shifts=len(shifts)
    )


def calculate_annotation_statistics(
    matched_ions: List[MatchedIon],
    theoretical_ions: List[TheoreticalIon],
    exp_mz: np.ndarray,
    exp_intensity: np.ndarray,
    peptide_length: int,
    min_sn: float = 5.0,
    use_resolved: bool = False,
) -> Dict:
    """
    Calculate comprehensive annotation statistics.

    Returns statistics similar to those reported by Annotator:
    - Sequence coverage (backbone bonds fragmented)
    - Peaks annotated (fraction of experimental peaks)
    - Intensity annotated (fraction of total intensity)
    - Fragments found (fraction of theoretical fragments)

    Args:
        matched_ions: List of matched ions
        theoretical_ions: List of all theoretical ions
        exp_mz: Experimental m/z array
        exp_intensity: Experimental intensity array
        peptide_length: Length of the peptide sequence
        min_sn: Minimum signal-to-noise ratio for bond coverage counting.
            Noise is estimated as the median peak intensity.
        use_resolved: If True, deduplicate so at most one ion per
            experimental peak contributes to the statistics (closest-ppm
            wins). Prevents a single peak from being counted multiple
            times when several theoretical ions happen to land on it
            (charge-reduced precursor envelope is a common case). Default
            False to preserve backward-compatible numbers.

    Returns:
        Dictionary with annotation statistics
    """
    # Optionally resolve to one ion per peak (closest-ppm wins)
    if use_resolved and matched_ions:
        _resolved = {}
        for ion in matched_ions:
            prev = _resolved.get(ion.exp_mz)
            if prev is None or abs(ion.mass_error_ppm) < abs(prev.mass_error_ppm):
                _resolved[ion.exp_mz] = ion
        matched_ions = list(_resolved.values())

    # Calculate sequence coverage (unique backbone positions fragmented)
    total_bonds = peptide_length - 1
    n_term_positions = set()  # b, c ions
    c_term_positions = set()  # y, z ions

    # Apply S/N filter for bond coverage to prevent noise peaks from
    # inflating coverage statistics
    noise_level = np.median(exp_intensity) if len(exp_intensity) > 0 else 1.0
    min_intensity = noise_level * min_sn

    for ion in matched_ions:
        if ion.exp_intensity < min_intensity:
            continue
        if ion.ion_type in ['b', 'c'] and ion.ion_number > 0:
            n_term_positions.add(ion.ion_number)
        elif ion.ion_type in ['y', 'z'] and ion.ion_number > 0:
            c_term_positions.add(ion.ion_number)

    # A bond is covered if we have either N-term or C-term fragment
    covered_bonds = set()
    for pos in n_term_positions:
        covered_bonds.add(pos)  # b_n/c_n covers bond after residue n
    for pos in c_term_positions:
        covered_bonds.add(peptide_length - pos)  # y_n/z_n covers bond before last n residues

    sequence_coverage = len(covered_bonds) / total_bonds if total_bonds > 0 else 0.0

    # Calculate peaks annotated
    matched_mz_set = set(ion.exp_mz for ion in matched_ions)
    peaks_annotated = len(matched_mz_set) / len(exp_mz) if len(exp_mz) > 0 else 0.0

    # Calculate intensity annotated
    total_intensity = np.sum(exp_intensity)
    matched_intensity = sum(ion.exp_intensity for ion in matched_ions)
    intensity_annotated = matched_intensity / total_intensity if total_intensity > 0 else 0.0

    # Calculate fragments found (theoretical fragments that were matched)
    # Count unique theoretical fragments (by type, number, charge)
    theo_backbone = [(ion.ion_type, ion.ion_number, ion.charge)
                     for ion in theoretical_ions
                     if ion.ion_type in ['b', 'y', 'c', 'z'] and not ion.neutral_loss]
    theo_backbone_unique = set(theo_backbone)

    matched_backbone = [(ion.ion_type, ion.ion_number, ion.charge)
                        for ion in matched_ions
                        if ion.ion_type in ['b', 'y', 'c', 'z'] and not ion.neutral_loss]
    matched_backbone_unique = set(matched_backbone)

    fragments_found = len(matched_backbone_unique) / len(theo_backbone_unique) if len(theo_backbone_unique) > 0 else 0.0

    return {
        'sequence_coverage': sequence_coverage,
        'sequence_coverage_bonds': f"{len(covered_bonds)}/{total_bonds}",
        'peaks_annotated': peaks_annotated,
        'peaks_annotated_count': f"{len(matched_mz_set)}/{len(exp_mz)}",
        'intensity_annotated': intensity_annotated,
        'fragments_found': fragments_found,
        'fragments_found_count': f"{len(matched_backbone_unique)}/{len(theo_backbone_unique)}",
        'n_term_coverage': n_term_positions,
        'c_term_coverage': c_term_positions,
    }


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_modifications_from_string(mod_string: str) -> List[Dict]:
    """
    Parse modification string from FragPipe format.

    Example: "N-term(229.1629),4S(528.2859),19K(229.1629)"
    """
    if not mod_string or mod_string == 'nan':
        return []

    mods = []
    for part in mod_string.split(','):
        part = part.strip()
        if not part or '(' not in part:
            continue

        try:
            mass_str = part[part.find('(')+1:part.find(')')]
            mass = float(mass_str)
            position_part = part[:part.find('(')]

            if position_part == 'N-term':
                mods.append({'position': 0, 'residue': 'N-term', 'mass': mass})
            elif position_part == 'C-term':
                mods.append({'position': -1, 'residue': 'C-term', 'mass': mass})
            else:
                # Parse position and residue
                pos = ''
                res = ''
                for char in position_part:
                    if char.isdigit():
                        pos += char
                    else:
                        res += char
                if pos and res:
                    mods.append({'position': int(pos), 'residue': res, 'mass': mass})
        except:
            continue

    return mods


# =============================================================================
# TEST / DEMO
# =============================================================================

if __name__ == "__main__":
    # Test with the example peptide from your data
    peptide = "AGYSQGATQYTQAQQTR"
    mod_string = "N-term(229.1629),4S(528.2859)"
    precursor_charge = 3

    print("=" * 70)
    print("Fragment Ion Calculator Test")
    print("=" * 70)
    print(f"Peptide: {peptide}")
    print(f"Modifications: {mod_string}")
    print(f"Precursor charge: +{precursor_charge}")

    # Parse modifications
    modifications = parse_modifications_from_string(mod_string)
    print(f"\nParsed modifications:")
    for mod in modifications:
        print(f"  Position {mod['position']}: {mod['residue']} +{mod['mass']:.4f}")

    # Create calculator
    calc = FragmentCalculator(peptide, modifications, precursor_charge)

    print(f"\nPrecursor mass: {calc.precursor_mass:.4f} Da")
    print(f"Precursor m/z: {calc.precursor_mz:.4f} (z={precursor_charge})")

    # Calculate all ions
    all_ions = calc.calculate_all_ions()

    print("\n" + "=" * 70)
    print("Theoretical Ions Summary")
    print("=" * 70)

    for ion_type, ions in all_ions.items():
        if ions:
            print(f"\n{ion_type} ions: {len(ions)}")
            # Show first few
            for ion in ions[:5]:
                print(f"  {ion.annotation}: m/z = {ion.mz:.4f}")
            if len(ions) > 5:
                print(f"  ... and {len(ions) - 5} more")

    # Test with example spectrum data
    print("\n" + "=" * 70)
    print("Peak Matching Test")
    print("=" * 70)

    # Load a test spectrum
    import pandas as pd
    spec_file = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/point_to_point_response/extracted_spectra_EThcD_HEK293T/Q96KR1_S195_18782.csv"

    try:
        spec_df = pd.read_csv(spec_file)
        exp_mz = spec_df['mz'].values
        exp_intensity = spec_df['intensity'].values

        # Get all theoretical ions
        all_theo_ions = calc.get_all_ions_flat()
        print(f"\nTotal theoretical ions: {len(all_theo_ions)}")

        # Match peaks
        matched = match_peaks(all_theo_ions, exp_mz, exp_intensity, tolerance_ppm=20.0)
        print(f"Matched ions: {len(matched)}")

        # Summarize by ion type
        from collections import Counter
        type_counts = Counter(ion.ion_type for ion in matched)
        print("\nMatched ions by type:")
        for ion_type, count in sorted(type_counts.items()):
            print(f"  {ion_type}: {count}")

        # Show top matched ions by intensity
        print("\nTop 10 matched ions by intensity:")
        matched_sorted = sorted(matched, key=lambda x: x.exp_intensity, reverse=True)
        for i, ion in enumerate(matched_sorted[:10], 1):
            print(f"  {i}. {ion.annotation}: theo={ion.mz:.4f}, exp={ion.exp_mz:.4f}, "
                  f"err={ion.mass_error_ppm:.1f}ppm, int={ion.exp_intensity:.0f}")

        # =================================================================
        # FALSE MATCH RATE CALCULATION
        # =================================================================
        print("\n" + "=" * 70)
        print("False Match Rate Calculation")
        print("(Spectrum shifting method from Schulte et al., Anal. Chem. 2025)")
        print("=" * 70)

        fmr = calculate_false_match_rate(
            all_theo_ions, exp_mz, exp_intensity,
            tolerance_ppm=20.0, shift_range=25.0, shift_step=1.0
        )

        print(f"\nTrue spectrum matches:")
        print(f"  Matched peaks: {fmr.matched_peaks}")
        print(f"  Matched intensity: {fmr.matched_intensity:.0f}")

        print(f"\nShifted spectrum matches (averaged over {fmr.n_shifts} shifts):")
        print(f"  Avg random peaks: {fmr.avg_random_peaks:.2f}")
        print(f"  Avg random intensity: {fmr.avg_random_intensity:.0f}")

        print(f"\nFalse Match Rates:")
        print(f"  FMR (peaks): {fmr.fmr_peaks*100:.1f}%")
        print(f"  FMR (intensity): {fmr.fmr_intensity*100:.1f}%")

        print("\nInterpretation:")
        print(f"  ~{fmr.fmr_peaks*100:.1f}% of matched peaks may be spurious")
        print(f"  ~{fmr.fmr_intensity*100:.1f}% of matched intensity may be spurious")

        # =================================================================
        # ANNOTATION STATISTICS
        # =================================================================
        print("\n" + "=" * 70)
        print("Annotation Statistics")
        print("=" * 70)

        stats = calculate_annotation_statistics(
            matched, all_theo_ions, exp_mz, exp_intensity, len(peptide)
        )

        print(f"\nSequence coverage: {stats['sequence_coverage_bonds']} "
              f"({stats['sequence_coverage']*100:.1f}%)")
        print(f"Peaks annotated: {stats['peaks_annotated_count']} "
              f"({stats['peaks_annotated']*100:.1f}%)")
        print(f"Intensity annotated: {stats['intensity_annotated']*100:.1f}%")
        print(f"Fragments found: {stats['fragments_found_count']} "
              f"({stats['fragments_found']*100:.1f}%)")

    except FileNotFoundError:
        print("Test spectrum file not found. Skipping peak matching test.")

    print("\n" + "=" * 70)
    print("Test complete!")
    print("=" * 70)
