#!/usr/bin/env python3
"""
Glycan Library for GlycoSpectrumAnnotator

This module contains glycan definitions, monosaccharide masses,
oxonium ions, and Y ion series calculations for various glycan types.

Supports:
- O-GlcNAc / O-GalNAc (simple O-glycans)
- N-glycans (high-mannose, complex, hybrid)
- Core structures and branching

Author: Longping Fu
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import re

# =============================================================================
# MONOSACCHARIDE MASSES (monoisotopic)
# =============================================================================

MONOSACCHARIDE_MASSES = {
    # Hexoses
    'Hex': 162.0528,        # Mannose, Glucose, Galactose
    'Man': 162.0528,        # Mannose (alias)
    'Glc': 162.0528,        # Glucose (alias)
    'Gal': 162.0528,        # Galactose (alias)

    # N-Acetylhexosamines
    'HexNAc': 203.0794,     # GlcNAc or GalNAc
    'GlcNAc': 203.0794,     # N-Acetylglucosamine (alias)
    'GalNAc': 203.0794,     # N-Acetylgalactosamine (alias)

    # Deoxyhexoses
    'Fuc': 146.0579,        # Fucose
    'dHex': 146.0579,       # Deoxyhexose (alias)

    # Sialic acids
    'NeuAc': 291.0954,      # N-Acetylneuraminic acid
    'Sia': 291.0954,        # Sialic acid (alias)
    'NeuGc': 307.0903,      # N-Glycolylneuraminic acid

    # Others
    'Xyl': 132.0423,        # Xylose
    'Pent': 132.0423,       # Pentose (alias)
    'HexA': 176.0321,       # Hexuronic acid (GlcA, IdoA)
    'GlcA': 176.0321,       # Glucuronic acid (alias)
    'Kdn': 250.0689,        # 2-keto-3-deoxy-nonulosonic acid
    'Sulfate': 79.9568,     # Sulfate modification
    'Phosphate': 79.9663,   # Phosphate modification
}

# =============================================================================
# GLYCAN TYPES AND COMPOSITIONS
# =============================================================================

@dataclass
class GlycanComposition:
    """Represents a glycan composition."""
    name: str
    composition: Dict[str, int]  # monosaccharide -> count
    mass: float
    glycan_type: str  # 'O-GlcNAc', 'O-GalNAc', 'N-glycan', 'O-glycan'

    @classmethod
    def from_string(cls, comp_string: str, name: str = "", glycan_type: str = "unknown"):
        """
        Parse glycan composition from string.

        Examples:
            "HexNAc1" -> {HexNAc: 1}
            "HexNAc2Hex3" -> {HexNAc: 2, Hex: 3}
            "HexNAc4Hex5Fuc1NeuAc2" -> {HexNAc: 4, Hex: 5, Fuc: 1, NeuAc: 2}
        """
        composition = {}
        pattern = r'([A-Za-z]+)(\d+)'

        for match in re.finditer(pattern, comp_string):
            mono = match.group(1)
            count = int(match.group(2))
            if mono in MONOSACCHARIDE_MASSES:
                composition[mono] = count

        mass = sum(MONOSACCHARIDE_MASSES[mono] * count
                   for mono, count in composition.items())

        return cls(
            name=name or comp_string,
            composition=composition,
            mass=mass,
            glycan_type=glycan_type
        )


# =============================================================================
# COMMON HUMAN O-GLYCAN COMPOSITIONS
# Based on MSFragger Human O-glycan database (12 common structures)
# Core structures: Core 1-4
# =============================================================================

O_GLYCAN_COMPOSITIONS = {
    # -------------------------------------------------------------------------
    # Simple O-GlcNAc/O-GalNAc (Tn antigen)
    # -------------------------------------------------------------------------
    'O-GlcNAc': GlycanComposition('O-GlcNAc', {'HexNAc': 1}, 203.0794, 'O-GlcNAc'),
    'O-GalNAc': GlycanComposition('O-GalNAc', {'HexNAc': 1}, 203.0794, 'O-GalNAc'),
    'Tn': GlycanComposition('Tn', {'HexNAc': 1}, 203.0794, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Core 1 structures (T antigen): Galβ1-3GalNAc
    # -------------------------------------------------------------------------
    'Core1': GlycanComposition('Core1', {'HexNAc': 1, 'Hex': 1}, 365.1322, 'O-glycan'),
    'T-antigen': GlycanComposition('T-antigen', {'HexNAc': 1, 'Hex': 1}, 365.1322, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Core 2 structures: Galβ1-3(GlcNAcβ1-6)GalNAc
    # -------------------------------------------------------------------------
    'Core2': GlycanComposition('Core2', {'HexNAc': 2, 'Hex': 1}, 568.2116, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Core 3 structures: GlcNAcβ1-3GalNAc
    # -------------------------------------------------------------------------
    'Core3': GlycanComposition('Core3', {'HexNAc': 2}, 406.1588, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Core 4 structures: GlcNAcβ1-3(GlcNAcβ1-6)GalNAc
    # -------------------------------------------------------------------------
    'Core4': GlycanComposition('Core4', {'HexNAc': 3}, 609.2382, 'O-glycan'),

    # -------------------------------------------------------------------------
    # MSFragger Human O-glycan database (12 common structures)
    # Format: HexNAc(n)Hex(m)Fuc(f)NeuAc(s)
    # -------------------------------------------------------------------------

    # 1. HexNAc(1) - Tn antigen (already defined above)

    # 2. HexNAc(1)Hex(1) - T antigen / Core 1 (already defined above)

    # 3. HexNAc(1)NeuAc(1) - Sialyl-Tn
    'HexNAc1NeuAc1': GlycanComposition('Sialyl-Tn', {'HexNAc': 1, 'NeuAc': 1}, 494.1748, 'O-glycan'),
    'Sialyl-Tn': GlycanComposition('Sialyl-Tn', {'HexNAc': 1, 'NeuAc': 1}, 494.1748, 'O-glycan'),

    # 4. HexNAc(2)Hex(1) - Core 2 (already defined above)

    # 5. HexNAc(1)Hex(1)NeuAc(1) - Sialyl-T / Sialyl-Core1
    'HexNAc1Hex1NeuAc1': GlycanComposition('Sialyl-T', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 1}, 656.2276, 'O-glycan'),
    'Sialyl-T': GlycanComposition('Sialyl-T', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 1}, 656.2276, 'O-glycan'),
    'Sialyl-Core1': GlycanComposition('Sialyl-Core1', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 1}, 656.2276, 'O-glycan'),

    # 6. HexNAc(2)Hex(2) - Extended Core 1 or Core 2
    'HexNAc2Hex2': GlycanComposition('HexNAc2Hex2', {'HexNAc': 2, 'Hex': 2}, 730.2644, 'O-glycan'),

    # 7. HexNAc(2)Hex(1)NeuAc(1) - Sialylated Core 2
    'HexNAc2Hex1NeuAc1': GlycanComposition('HexNAc2Hex1NeuAc1', {'HexNAc': 2, 'Hex': 1, 'NeuAc': 1}, 859.3070, 'O-glycan'),

    # 8. HexNAc(1)Hex(1)NeuAc(2) - Disialyl-T
    'HexNAc1Hex1NeuAc2': GlycanComposition('Disialyl-T', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 2}, 947.3230, 'O-glycan'),
    'Disialyl-T': GlycanComposition('Disialyl-T', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 2}, 947.3230, 'O-glycan'),

    # 9. HexNAc(2)Hex(2)NeuAc(1) - Monosialylated extended
    'HexNAc2Hex2NeuAc1': GlycanComposition('HexNAc2Hex2NeuAc1', {'HexNAc': 2, 'Hex': 2, 'NeuAc': 1}, 1021.3598, 'O-glycan'),

    # 10. HexNAc(2)Hex(2)Fuc(1)NeuAc(1) - Fucosylated monosialylated
    'HexNAc2Hex2Fuc1NeuAc1': GlycanComposition('HexNAc2Hex2Fuc1NeuAc1', {'HexNAc': 2, 'Hex': 2, 'Fuc': 1, 'NeuAc': 1}, 1167.4177, 'O-glycan'),

    # 11. HexNAc(2)Hex(2)NeuAc(2) - Disialylated extended
    'HexNAc2Hex2NeuAc2': GlycanComposition('HexNAc2Hex2NeuAc2', {'HexNAc': 2, 'Hex': 2, 'NeuAc': 2}, 1312.4552, 'O-glycan'),

    # 12. HexNAc(2)Hex(2)Fuc(1)NeuAc(2) - Fucosylated disialylated
    'HexNAc2Hex2Fuc1NeuAc2': GlycanComposition('HexNAc2Hex2Fuc1NeuAc2', {'HexNAc': 2, 'Hex': 2, 'Fuc': 1, 'NeuAc': 2}, 1458.5131, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Additional common O-glycans with NeuGc (non-human, but for completeness)
    # -------------------------------------------------------------------------
    'HexNAc1NeuGc1': GlycanComposition('HexNAc1NeuGc1', {'HexNAc': 1, 'NeuGc': 1}, 510.1697, 'O-glycan'),
    'HexNAc1Hex1NeuGc1': GlycanComposition('HexNAc1Hex1NeuGc1', {'HexNAc': 1, 'Hex': 1, 'NeuGc': 1}, 672.2225, 'O-glycan'),

    # -------------------------------------------------------------------------
    # Phosphorylated and Sulfated O-glycans (common in your lab's data)
    # -------------------------------------------------------------------------
    'HexNAc1-Phosphate': GlycanComposition('HexNAc1-Phosphate', {'HexNAc': 1, 'Phosphate': 1}, 283.0457, 'O-glycan'),
    'HexNAc1Hex1-Phosphate': GlycanComposition('HexNAc1Hex1-Phosphate', {'HexNAc': 1, 'Hex': 1, 'Phosphate': 1}, 445.0985, 'O-glycan'),
    'HexNAc1-Sulfate': GlycanComposition('HexNAc1-Sulfate', {'HexNAc': 1, 'Sulfate': 1}, 283.0362, 'O-glycan'),
    'HexNAc1Hex1-Sulfate': GlycanComposition('HexNAc1Hex1-Sulfate', {'HexNAc': 1, 'Hex': 1, 'Sulfate': 1}, 445.0890, 'O-glycan'),
    'HexNAc2Hex2-Sulfate': GlycanComposition('HexNAc2Hex2-Sulfate', {'HexNAc': 2, 'Hex': 2, 'Sulfate': 1}, 810.2212, 'O-glycan'),
}

# Common N-glycan compositions
N_GLYCAN_COMPOSITIONS = {
    # Trimannosyl core (all N-glycans share this)
    'Core': GlycanComposition('N-glycan Core', {'HexNAc': 2, 'Hex': 3}, 892.3170, 'N-glycan'),

    # High-mannose series
    'Man5': GlycanComposition('Man5', {'HexNAc': 2, 'Hex': 5}, 1216.4226, 'N-glycan'),
    'Man6': GlycanComposition('Man6', {'HexNAc': 2, 'Hex': 6}, 1378.4754, 'N-glycan'),
    'Man7': GlycanComposition('Man7', {'HexNAc': 2, 'Hex': 7}, 1540.5282, 'N-glycan'),
    'Man8': GlycanComposition('Man8', {'HexNAc': 2, 'Hex': 8}, 1702.5810, 'N-glycan'),
    'Man9': GlycanComposition('Man9', {'HexNAc': 2, 'Hex': 9}, 1864.6338, 'N-glycan'),

    # Complex N-glycans (biantennary)
    'A2': GlycanComposition('A2', {'HexNAc': 4, 'Hex': 5}, 1622.5814, 'N-glycan'),
    'A2F': GlycanComposition('A2F', {'HexNAc': 4, 'Hex': 5, 'Fuc': 1}, 1768.6393, 'N-glycan'),
    'A2G2': GlycanComposition('A2G2', {'HexNAc': 4, 'Hex': 5}, 1622.5814, 'N-glycan'),
    'A2G2F': GlycanComposition('A2G2F', {'HexNAc': 4, 'Hex': 5, 'Fuc': 1}, 1768.6393, 'N-glycan'),
    'A2G2S1': GlycanComposition('A2G2S1', {'HexNAc': 4, 'Hex': 5, 'NeuAc': 1}, 1913.6768, 'N-glycan'),
    'A2G2S2': GlycanComposition('A2G2S2', {'HexNAc': 4, 'Hex': 5, 'NeuAc': 2}, 2204.7722, 'N-glycan'),
    'A2G2FS1': GlycanComposition('A2G2FS1', {'HexNAc': 4, 'Hex': 5, 'Fuc': 1, 'NeuAc': 1}, 2059.7347, 'N-glycan'),
    'A2G2FS2': GlycanComposition('A2G2FS2', {'HexNAc': 4, 'Hex': 5, 'Fuc': 1, 'NeuAc': 2}, 2350.8301, 'N-glycan'),
}

# =============================================================================
# OXONIUM IONS (Glycan Diagnostic B-ions)
# =============================================================================

# Extended oxonium ion library
OXONIUM_IONS_EXTENDED = {
    # HexNAc (GlcNAc/GalNAc) series
    'HexNAc': 204.0867,              # [HexNAc + H]+
    'HexNAc-H2O': 186.0761,          # [HexNAc - H2O + H]+
    'HexNAc-2H2O': 168.0655,         # [HexNAc - 2H2O + H]+
    'HexNAc-CH2O': 174.0761,         # Cross-ring fragment
    'HexNAc-C2H4O2': 144.0655,       # Cross-ring fragment
    'HexNAc-C2H6O3': 126.0550,       # Cross-ring fragment
    'HexNAc_138': 138.0550,          # Common fragment

    # Hex (Man/Glc/Gal) series
    'Hex': 163.0601,                 # [Hex + H]+
    'Hex-H2O': 145.0495,             # [Hex - H2O + H]+
    'Hex-2H2O': 127.0390,            # [Hex - 2H2O + H]+

    # NeuAc (Sialic acid) series
    'NeuAc': 292.1027,               # [NeuAc + H]+
    'NeuAc-H2O': 274.0921,           # [NeuAc - H2O + H]+
    'NeuAc-2H2O': 256.0816,          # [NeuAc - 2H2O + H]+
    'NeuAc-CO2': 248.1129,           # [NeuAc - CO2 + H]+
    'NeuAc-CH4O2': 260.0765,         # Cross-ring fragment

    # Fuc (Fucose) series
    'Fuc': 147.0652,                 # [Fuc + H]+
    'Fuc-H2O': 129.0546,             # [Fuc - H2O + H]+

    # Disaccharide oxonium ions
    'HexNAc-Hex': 366.1395,          # [HexNAc + Hex + H]+
    'HexNAc-Hex-H2O': 348.1289,      # [HexNAc + Hex - H2O + H]+
    'HexNAc-HexNAc': 407.1660,       # [2 HexNAc + H]+
    'Hex-Hex': 325.1129,             # [2 Hex + H]+
    'HexNAc-NeuAc': 495.1821,        # [HexNAc + NeuAc + H]+
    'Hex-NeuAc': 454.1556,           # [Hex + NeuAc + H]+

    # Trisaccharide oxonium ions
    'HexNAc-Hex-NeuAc': 657.2349,    # Common sialylated fragment
    'HexNAc-2Hex': 528.1923,         # HexNAc + 2Hex (high-mannose arm)

    # Tetrasaccharide+ oxonium ions (high-mannose glycan arm fragments)
    'HexNAc-3Hex': 690.2451,         # HexNAc + 3Hex
    'HexNAc-4Hex': 852.2979,         # HexNAc + 4Hex

    # TMT-modified oxonium ions (for TMT-labeled glycopeptides)
    'HexNAc-TMT': 529.2937,          # [HexNAc + TMT + H]+
    'HexNAc-H2O-TMT': 511.2831,      # [HexNAc - H2O + TMT + H]+
}

# Subset for O-GlcNAc/O-GalNAc analysis
OXONIUM_IONS_O_GLCNAC = {
    'HexNAc': 204.0867,
    'HexNAc-H2O': 186.0761,
    'HexNAc-2H2O': 168.0655,
    'HexNAc_138': 138.0550,
    'HexNAc_126': 126.0550,
    'HexNAc-TMT': 529.2937,
}

# Subset for N-glycan analysis
OXONIUM_IONS_N_GLYCAN = {
    'HexNAc': 204.0867,
    'HexNAc-H2O': 186.0761,
    'HexNAc-2H2O': 168.0655,
    'HexNAc_138': 138.0550,
    'Hex': 163.0601,
    'Hex-H2O': 145.0495,
    'NeuAc': 292.1027,
    'NeuAc-H2O': 274.0921,
    'Fuc': 147.0652,
    'HexNAc-Hex': 366.1395,
    'HexNAc-Hex-H2O': 348.1289,
    'HexNAc-NeuAc': 495.1821,
    'HexNAc-2Hex': 528.1923,         # HexNAc + 2Hex (high-mannose arm)
    'HexNAc-Hex-NeuAc': 657.2349,
    'HexNAc-3Hex': 690.2451,         # HexNAc + 3Hex (high-mannose arm)
    'HexNAc-4Hex': 852.2979,         # HexNAc + 4Hex (high-mannose arm)
}

# Extended set for O-glycan analysis (mucin-type cores 1/2/3/4 + sialylation/fucosylation).
# Excludes high-mannose-style ions like HexNAc-2Hex / HexNAc-3Hex / HexNAc-4Hex which
# do not occur in typical mucin-type O-glycan fragmentation.
OXONIUM_IONS_O_GLYCAN = {
    'HexNAc_126':       126.0550,
    'HexNAc_138':       138.0550,
    'HexNAc_144':       144.0656,
    'Fuc':              147.0652,
    'Hex':              163.0601,
    'HexNAc-2H2O':      168.0655,
    'HexNAc-H2O':       186.0761,
    'HexNAc':           204.0867,
    'NeuAc-H2O':        274.0921,
    'NeuAc':            292.1027,
    'HexNAc+':          300.1308,   # HexNAc + H2O + H
    'HexNAc-Fuc':       350.1340,   # core 1 + Fuc fragment
    'HexNAc-Hex':       366.1395,   # T-antigen
    'HexNAc-Hex-Fuc':   512.1974,
    'HexNAc-NeuAc':     495.1821,   # sialyl-Tn
    'HexNAc-Hex-NeuAc': 657.2349,   # sialyl-T
    'HexNAc-TMT':       529.2937,
}

# =============================================================================
# Y ION SERIES GENERATION
# =============================================================================

def generate_y_ion_series(
    glycan: GlycanComposition,
    peptide_mass: float,
    include_water_loss: bool = True,
    tmt_on_glycan: bool = False,
    tmt_mass: float = 229.1629
) -> Dict[str, float]:
    """
    Generate Y ion series for a glycan composition.

    Y ions represent the peptide backbone with varying amounts of glycan attached:
    - Y0: Peptide only (complete glycan loss)
    - Y1, Y2, ...: Peptide + increasing glycan fragments
    - Y(intact): Full glycopeptide

    Args:
        glycan: GlycanComposition object
        peptide_mass: Mass of the peptide backbone (without glycan)
        include_water_loss: Include Y-H2O variants
        tmt_on_glycan: Whether TMT tag is on the glycan
        tmt_mass: Mass of TMT tag

    Returns:
        Dictionary of Y ion names to neutral masses
    """
    y_ions = {}
    glycan_mass = glycan.mass

    if tmt_on_glycan:
        glycan_mass += tmt_mass

    # Y0 - peptide only (complete glycan loss)
    y_ions['Y0'] = peptide_mass
    if include_water_loss:
        y_ions['Y0-H2O'] = peptide_mass - 18.0106

    # For simple O-glycans (single monosaccharide)
    if glycan.glycan_type in ['O-GlcNAc', 'O-GalNAc']:
        # Y1 = intact glycopeptide
        y_ions['Y1'] = peptide_mass + glycan_mass
        if include_water_loss:
            y_ions['Y1-H2O'] = peptide_mass + glycan_mass - 18.0106

    # For complex O-glycans and N-glycans
    else:
        comp = glycan.composition

        # Build Y ion ladder based on composition
        # Strategy: remove monosaccharides one at a time from non-reducing end

        # Y ions by number of monosaccharides attached
        total_mono = sum(comp.values())

        for i in range(1, total_mono + 1):
            # Generate representative Y ions
            # This is a simplified model - full implementation would track topology

            if i == 1 and 'HexNAc' in comp:
                # Y1 = peptide + 1 HexNAc (reducing end)
                y_ions[f'Y1'] = peptide_mass + MONOSACCHARIDE_MASSES['HexNAc']

            if i == 2 and comp.get('HexNAc', 0) >= 2:
                # Y2 = peptide + 2 HexNAc (chitobiose core for N-glycans)
                y_ions[f'Y2'] = peptide_mass + 2 * MONOSACCHARIDE_MASSES['HexNAc']

            # For N-glycans, add core structure Y ions
            if glycan.glycan_type == 'N-glycan':
                # Trimannosyl core
                if comp.get('HexNAc', 0) >= 2 and comp.get('Hex', 0) >= 3:
                    core_mass = 2 * MONOSACCHARIDE_MASSES['HexNAc'] + 3 * MONOSACCHARIDE_MASSES['Hex']
                    y_ions['Y(core)'] = peptide_mass + core_mass

                    # Add fucose to core if present
                    if comp.get('Fuc', 0) >= 1:
                        y_ions['Y(core+F)'] = peptide_mass + core_mass + MONOSACCHARIDE_MASSES['Fuc']

        # Y(intact) = full glycopeptide
        y_ions['Y(intact)'] = peptide_mass + glycan_mass

        # Add water loss variants for major Y ions
        if include_water_loss:
            for name, mass in list(y_ions.items()):
                if not name.endswith('-H2O'):
                    y_ions[f'{name}-H2O'] = mass - 18.0106

    return y_ions


def generate_n_glycan_y_ions(
    glycan_composition: Dict[str, int],
    peptide_mass: float,
    include_fucose_variants: bool = True
) -> Dict[str, float]:
    """
    Generate comprehensive Y ion series for N-glycans.

    N-glycan Y ions follow the fragmentation pattern:
    - Loss from non-reducing end (antenna)
    - Core structure usually retained

    Args:
        glycan_composition: Dict of monosaccharide counts
        peptide_mass: Mass of peptide backbone
        include_fucose_variants: Include +/- Fuc variants

    Returns:
        Dictionary of Y ion names to neutral masses
    """
    y_ions = {}

    n_hexnac = glycan_composition.get('HexNAc', 0)
    n_hex = glycan_composition.get('Hex', 0)
    n_fuc = glycan_composition.get('Fuc', 0)
    n_neuac = glycan_composition.get('NeuAc', 0)

    hex_mass = MONOSACCHARIDE_MASSES['Hex']
    hexnac_mass = MONOSACCHARIDE_MASSES['HexNAc']
    fuc_mass = MONOSACCHARIDE_MASSES['Fuc']
    neuac_mass = MONOSACCHARIDE_MASSES['NeuAc']

    # Y0 - peptide only
    y_ions['Y0'] = peptide_mass

    # Y1 - peptide + 1 HexNAc (reducing end GlcNAc)
    if n_hexnac >= 1:
        y_ions['Y1'] = peptide_mass + hexnac_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y1F'] = peptide_mass + hexnac_mass + fuc_mass

    # Y2 - peptide + chitobiose (2 HexNAc)
    if n_hexnac >= 2:
        y_ions['Y2'] = peptide_mass + 2 * hexnac_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y2F'] = peptide_mass + 2 * hexnac_mass + fuc_mass

    # Y3 - peptide + chitobiose + 1 Man
    if n_hexnac >= 2 and n_hex >= 1:
        y_ions['Y3'] = peptide_mass + 2 * hexnac_mass + hex_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y3F'] = peptide_mass + 2 * hexnac_mass + hex_mass + fuc_mass

    # Y4 - peptide + chitobiose + 2 Man
    if n_hexnac >= 2 and n_hex >= 2:
        y_ions['Y4'] = peptide_mass + 2 * hexnac_mass + 2 * hex_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y4F'] = peptide_mass + 2 * hexnac_mass + 2 * hex_mass + fuc_mass

    # Y(core) - trimannosyl core (HexNAc2Hex3)
    if n_hexnac >= 2 and n_hex >= 3:
        core_mass = 2 * hexnac_mass + 3 * hex_mass
        y_ions['Y(core)'] = peptide_mass + core_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y(core)F'] = peptide_mass + core_mass + fuc_mass

    # For complex N-glycans, add antenna Y ions
    if n_hexnac >= 3 and n_hex >= 3:
        # One antenna
        y_ions['Y(core+1arm)'] = peptide_mass + 3 * hexnac_mass + 4 * hex_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y(core+1arm)F'] = peptide_mass + 3 * hexnac_mass + 4 * hex_mass + fuc_mass

    if n_hexnac >= 4 and n_hex >= 5:
        # Both antennae (biantennary)
        y_ions['Y(bi)'] = peptide_mass + 4 * hexnac_mass + 5 * hex_mass
        if n_fuc >= 1 and include_fucose_variants:
            y_ions['Y(bi)F'] = peptide_mass + 4 * hexnac_mass + 5 * hex_mass + fuc_mass

    # Add sialylated variants
    if n_neuac >= 1:
        # Add NeuAc to major Y ions
        for name, mass in list(y_ions.items()):
            if 'core' in name or 'bi' in name:
                y_ions[f'{name}S1'] = mass + neuac_mass
        if n_neuac >= 2:
            for name, mass in list(y_ions.items()):
                if ('core' in name or 'bi' in name) and 'S' not in name:
                    y_ions[f'{name}S2'] = mass + 2 * neuac_mass

    # Full glycopeptide
    total_mass = (n_hexnac * hexnac_mass + n_hex * hex_mass +
                  n_fuc * fuc_mass + n_neuac * neuac_mass)
    y_ions['Y(intact)'] = peptide_mass + total_mass

    return y_ions


# =============================================================================
# GLYCAN IDENTIFICATION FROM MASS
# =============================================================================

def identify_glycan_from_mass(
    mass: float,
    tolerance_da: float = 0.1,
    glycan_type: str = 'all'
) -> List[Tuple[str, GlycanComposition, float]]:
    """
    Identify potential glycan compositions from a mass.

    Args:
        mass: Observed glycan mass
        tolerance_da: Mass tolerance in Daltons
        glycan_type: 'O-glycan', 'N-glycan', or 'all'

    Returns:
        List of (name, composition, mass_error) tuples
    """
    matches = []

    if glycan_type in ['O-glycan', 'all']:
        for name, glycan in O_GLYCAN_COMPOSITIONS.items():
            error = abs(mass - glycan.mass)
            if error <= tolerance_da:
                matches.append((name, glycan, error))

    if glycan_type in ['N-glycan', 'all']:
        for name, glycan in N_GLYCAN_COMPOSITIONS.items():
            error = abs(mass - glycan.mass)
            if error <= tolerance_da:
                matches.append((name, glycan, error))

    return sorted(matches, key=lambda x: x[2])


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_glycan_mass(composition_string: str) -> float:
    """
    Calculate mass from composition string.

    Example: "HexNAc2Hex5" -> 1216.4226
    """
    glycan = GlycanComposition.from_string(composition_string)
    return glycan.mass


def parse_proforma_glycan(proforma_string: str) -> Optional[GlycanComposition]:
    """
    Parse ProForma 2.0 glycan notation.

    Examples:
        "[Glycan:HexNAc1]" -> GlycanComposition for O-GlcNAc
        "[Glycan:HexNAc2Hex5]" -> GlycanComposition for Man5
    """
    match = re.search(r'\[Glycan:([A-Za-z0-9]+)\]', proforma_string)
    if match:
        comp_string = match.group(1)
        return GlycanComposition.from_string(comp_string)
    return None


# =============================================================================
# CROSSLINKER SUPPORT
# =============================================================================
# MS-cleavable crosslinkers (DSSO, DSBSO) and non-cleavable (BS3)
# Reference: Schulte et al. Anal. Chem. 2025, 97, 23120-23130
# =============================================================================

@dataclass
class Crosslinker:
    """
    Represents a chemical crosslinker for XL-MS experiments.

    Attributes:
        name: Crosslinker name
        formula: Chemical formula
        intact_mass: Mass of intact crosslinker (after reaction, losing NHS groups)
        spacer_length: Spacer arm length in Angstroms
        cleavable: Whether MS-cleavable
        reactive_groups: What residues it reacts with (e.g., 'NHS' for Lys)
        stub_masses: Dict of stub names to masses (for cleavable crosslinkers)
        diagnostic_ions: Dict of diagnostic ion names to m/z values
    """
    name: str
    formula: str
    intact_mass: float
    spacer_length: float
    cleavable: bool
    reactive_groups: str
    stub_masses: Dict[str, float]
    diagnostic_ions: Dict[str, float]


# =============================================================================
# CROSSLINKER DEFINITIONS
# =============================================================================

# DSSO - Disuccinimidyl sulfoxide (MS-cleavable, NHS-ester)
# Thermo Fisher Scientific, commonly used MS-cleavable crosslinker
# Cleaves at C-S bonds adjacent to sulfoxide, producing alkene and thiol/sulfenic stubs
DSSO = Crosslinker(
    name='DSSO',
    formula='C14H18N2O9S',
    intact_mass=158.0038,  # Mass added to crosslinked peptide pair (after losing 2x NHS)
    spacer_length=10.1,    # Angstroms
    cleavable=True,
    reactive_groups='NHS',  # Reacts with Lys, protein N-terminus
    stub_masses={
        'alkene': 54.0106,      # C3H2O - Alkene stub (A)
        'thiol': 85.9826,       # C3H2OS - Unsaturated thiol stub (T)
        'sulfenic': 103.9932,   # C3H4O2S - Sulfenic acid stub (S)
    },
    diagnostic_ions={
        # Mass difference between A and T stubs = 31.97 Da
        'A-T_diff': 31.9720,
        # Reporter ions in low mass region
        'reporter_104': 104.0170,
        'reporter_122': 122.0276,
    }
)

# DSBSO - Disuccinimidyl bis-sulfoxide (MS-cleavable, enrichable)
# Contains azide/alkyne for click chemistry enrichment + acid-cleavable site
# Similar sulfoxide cleavage mechanism to DSSO
# Reference: Kao et al. Mol Cell Proteomics 2011; Jiang et al. Angew Chem 2026
DSBSO = Crosslinker(
    name='DSBSO',
    formula='C11H16N2O6S2',
    intact_mass=308.0388,  # Full crosslinker mass (C11H16O6S2 portion)
    spacer_length=12.5,    # Angstroms (approximate)
    cleavable=True,
    reactive_groups='NHS',
    stub_masses={
        # DSBSO has similar sulfoxide cleavage chemistry
        # Stub masses may differ from DSSO due to different spacer
        'alkene': 54.0106,      # C3H2O - Alkene stub (A)
        'thiol': 85.9826,       # Unsaturated thiol stub (T)
        'sulfenic': 103.9932,   # Sulfenic acid stub (S)
        # Extended fragments specific to DSBSO
        'ETHMP': 226.0144,      # Extended thiol fragment
    },
    diagnostic_ions={
        'A-T_diff': 31.9720,
        # DSBSO-specific diagnostic ions
        'D12_ETHMP_alkene_diff': 226.01,  # Distance between D12 and ETHMP+alkene
    }
)

# BS3 - Bis(sulfosuccinimidyl)suberate (non-cleavable, NHS-ester)
# Classic non-cleavable crosslinker
BS3 = Crosslinker(
    name='BS3',
    formula='C16H18N2O14S2',
    intact_mass=138.0681,  # Suberate spacer mass (C8H12O2) after NHS loss
    spacer_length=11.4,    # Angstroms
    cleavable=False,
    reactive_groups='NHS',
    stub_masses={},  # Non-cleavable, no stubs
    diagnostic_ions={}
)

# DSS - Disuccinimidyl suberate (non-cleavable, membrane permeable version of BS3)
DSS = Crosslinker(
    name='DSS',
    formula='C16H20N2O8',
    intact_mass=138.0681,  # Same spacer as BS3
    spacer_length=11.4,
    cleavable=False,
    reactive_groups='NHS',
    stub_masses={},
    diagnostic_ions={}
)

# Dictionary of all crosslinkers
CROSSLINKERS = {
    'DSSO': DSSO,
    'DSBSO': DSBSO,
    'BS3': BS3,
    'DSS': DSS,
}


# =============================================================================
# CROSSLINKED PEPTIDE FRAGMENTATION
# =============================================================================

def generate_crosslink_fragments(
    peptide1_mass: float,
    peptide2_mass: float,
    crosslinker: Crosslinker,
    precursor_charge: int,
    include_neutral_losses: bool = True
) -> Dict[str, Dict[str, float]]:
    """
    Generate theoretical fragments for a crosslinked peptide pair.

    For MS-cleavable crosslinkers (DSSO, DSBSO):
    - Cleaves at sulfoxide C-S bonds
    - Produces signature doublet pairs (αA/βT and αT/βA)
    - Each peptide retains one stub (alkene, thiol, or sulfenic acid)

    Args:
        peptide1_mass: Neutral mass of first peptide
        peptide2_mass: Neutral mass of second peptide
        crosslinker: Crosslinker object
        precursor_charge: Charge state of precursor
        include_neutral_losses: Include H2O and NH3 losses

    Returns:
        Dictionary with 'peptide1', 'peptide2', and 'diagnostic' fragments
    """
    fragments = {
        'peptide1': {},
        'peptide2': {},
        'diagnostic': {},
        'precursor': {},
    }

    PROTON = 1.007276
    H2O = 18.0106
    NH3 = 17.0265

    if crosslinker.cleavable:
        # MS-cleavable crosslinker fragmentation
        alkene = crosslinker.stub_masses.get('alkene', 0)
        thiol = crosslinker.stub_masses.get('thiol', 0)
        sulfenic = crosslinker.stub_masses.get('sulfenic', 0)

        # Peptide 1 with different stubs
        fragments['peptide1']['α-A'] = peptide1_mass + alkene  # Alkene stub
        fragments['peptide1']['α-T'] = peptide1_mass + thiol   # Thiol stub
        fragments['peptide1']['α-S'] = peptide1_mass + sulfenic  # Sulfenic stub

        # Peptide 2 with different stubs
        fragments['peptide2']['β-A'] = peptide2_mass + alkene
        fragments['peptide2']['β-T'] = peptide2_mass + thiol
        fragments['peptide2']['β-S'] = peptide2_mass + sulfenic

        # Add charged versions
        for charge in range(1, min(precursor_charge, 4)):
            for name, mass in list(fragments['peptide1'].items()):
                fragments['peptide1'][f'{name}+{charge}'] = (mass + charge * PROTON) / charge
            for name, mass in list(fragments['peptide2'].items()):
                fragments['peptide2'][f'{name}+{charge}'] = (mass + charge * PROTON) / charge

        # Neutral losses
        if include_neutral_losses:
            for pep_key in ['peptide1', 'peptide2']:
                for name, mass in list(fragments[pep_key].items()):
                    if '+' not in name:  # Only for neutral masses
                        fragments[pep_key][f'{name}-H2O'] = mass - H2O
                        fragments[pep_key][f'{name}-NH3'] = mass - NH3

        # Diagnostic ions
        for ion_name, ion_mass in crosslinker.diagnostic_ions.items():
            fragments['diagnostic'][ion_name] = ion_mass

    else:
        # Non-cleavable crosslinker - intact crosslinked peptide
        intact_mass = peptide1_mass + peptide2_mass + crosslinker.intact_mass

        fragments['precursor']['intact'] = intact_mass

        # Charged precursors
        for charge in range(1, precursor_charge + 1):
            fragments['precursor'][f'intact+{charge}'] = (intact_mass + charge * PROTON) / charge

    return fragments


def identify_crosslink_stubs(
    observed_masses: List[float],
    peptide_mass: float,
    crosslinker: Crosslinker,
    tolerance_da: float = 0.02
) -> List[Tuple[str, float, float]]:
    """
    Identify crosslinker stub modifications on a peptide.

    Looks for signature mass shifts indicating crosslinker cleavage.

    Args:
        observed_masses: List of observed fragment masses
        peptide_mass: Expected peptide mass without crosslinker
        crosslinker: Crosslinker object
        tolerance_da: Mass tolerance

    Returns:
        List of (stub_name, observed_mass, mass_error) tuples
    """
    matches = []

    if not crosslinker.cleavable:
        return matches

    for stub_name, stub_mass in crosslinker.stub_masses.items():
        expected = peptide_mass + stub_mass

        for obs_mass in observed_masses:
            error = abs(obs_mass - expected)
            if error <= tolerance_da:
                matches.append((stub_name, obs_mass, error))

    return matches


def calculate_crosslink_fmr(
    matched_peaks: int,
    total_peaks: int,
    matched_intensity: float,
    total_intensity: float,
    shift_range: int = 50,
    pi_offset: float = 3.14159
) -> Dict[str, float]:
    """
    Calculate False Match Rate for crosslink spectrum annotation.

    Uses the spectrum shifting method from Schulte et al.:
    - Shift spectrum by π ± 25 Th in 1 Th steps
    - π offset prevents isotope pattern false matches
    - FMR = average(shifted matches) / unshifted matches

    Args:
        matched_peaks: Number of matched peaks in original spectrum
        total_peaks: Total number of peaks
        matched_intensity: Sum of matched peak intensities
        total_intensity: Sum of all peak intensities
        shift_range: Range of shifts (default ±25 Th = 50 shifts)
        pi_offset: Offset to avoid isotope patterns (default π)

    Returns:
        Dictionary with FMR estimates for peaks and intensity
    """
    # This is a placeholder for the actual FMR calculation
    # Real implementation requires the full spectrum data
    return {
        'fmr_peaks': matched_peaks / total_peaks if total_peaks > 0 else 0,
        'fmr_intensity': matched_intensity / total_intensity if total_intensity > 0 else 0,
        'note': 'Full FMR calculation requires spectrum shifting - implement in annotator'
    }


# =============================================================================
# PROFORMA CROSSLINK NOTATION
# =============================================================================

def parse_proforma_crosslink(proforma_string: str) -> Optional[Tuple[str, str, str]]:
    """
    Parse ProForma 2.0 crosslink notation.

    Examples:
        "KC[L-cystine#XL1]M//GC[#XL1]V" -> ('KC[...]M', 'GC[...]V', 'L-cystine')
        "PEPTIDEK[DSSO#XL1]//ANOTHERK[#XL1]" -> peptide pair with DSSO

    Returns:
        Tuple of (peptide1, peptide2, crosslinker_name) or None
    """
    # Check for chimeric peptide separator
    if '//' not in proforma_string:
        return None

    parts = proforma_string.split('//')
    if len(parts) != 2:
        return None

    peptide1, peptide2 = parts

    # Extract crosslinker name from first occurrence
    xl_match = re.search(r'\[([A-Za-z0-9-]+)#XL\d+\]', peptide1)
    if xl_match:
        crosslinker_name = xl_match.group(1)
        return (peptide1, peptide2, crosslinker_name)

    return None


# ============================================================================
# N-glycan HCD fragment generators (moved from nglyco_annotate_spectra 2026-07-06)
# Composition-aware remainder b/y masses + diagnostic B-ions, kept in the shared
# engine so every N-glyco annotator/test draws from one implementation (mirrors
# generate_n_glycan_y_ions above). Self-contained: local monosaccharide constants
# preserve exact numeric behavior.
# ============================================================================

def _classify_glycan_type_simple(glycan_comp):
    """Lightweight glycan type classification for fragment generation.

    Returns: 'high_mannose', 'paucimannose', 'complex', 'hybrid', or 'unknown'.
    """
    n = glycan_comp.get('HexNAc', 0)
    h = glycan_comp.get('Hex', 0)
    sia = glycan_comp.get('NeuAc', 0) + glycan_comp.get('NeuGc', 0)

    if n < 2:
        return 'unknown'
    if n == 2 and h >= 5 and sia == 0:
        return 'high_mannose'
    if n == 2 and h <= 4:
        return 'paucimannose'
    if n == 3 and h >= 5:
        return 'hybrid'
    if n >= 3:
        return 'complex'
    return 'unknown'


def generate_n_glycan_remainders(glycan_comp):
    """
    Generate biosynthetically valid partial glycan masses for N-glycan HCD.

    Applies glycan-type-aware fragmentation logic based on the N-glycan
    biosynthetic tree structure:

    High-mannose (N2H5-9):
        Sequential mannose loss from non-reducing end.  All intermediates
        N2H3→N2H4→...→N2H(n_hex) are valid (linear mannose tree).

    Complex/hybrid (N3+ H3+):
        Antennae are GlcNAc-Gal(-NeuAc) units connected to the trimannosyl
        core through GlcNAc bridges.  Gal CANNOT exist without its preceding
        GlcNAc (Riley JPR 2020; StrucGAP biosynthetic rules).  Valid
        intermediates above the core must contain complete or partially-stripped
        antenna units (GlcNAc, then GlcNAc-Gal, then GlcNAc-Gal-NeuAc).

    Returns dict: {label: mass}. Label '' = bare, 'intact' = full glycan.
    """
    N = 203.0794   # HexNAc
    H = 162.0528   # Hex (Man, Gal, Glc)
    F = 146.0579   # Fuc
    A = 291.0954   # NeuAc
    G = 307.0903   # NeuGc

    n_N = glycan_comp.get('HexNAc', 0)
    n_H = glycan_comp.get('Hex', 0)
    n_F = glycan_comp.get('Fuc', 0)
    n_A = glycan_comp.get('NeuAc', 0)
    n_G = glycan_comp.get('NeuGc', 0)

    total = n_N * N + n_H * H + n_F * F + n_A * A + n_G * G

    rem = {'': 0.0}  # Y0 / bare (complete glycan loss)

    gtype = _classify_glycan_type_simple(glycan_comp)

    # =================================================================
    # Chitobiose core ladder — universal for all N-glycans
    # GlcNAc(β1-4)GlcNAc—Asn
    # =================================================================
    if n_N >= 1:
        rem['+N'] = N                          # Y1: reducing-end GlcNAc
    if n_N >= 2:
        rem['+N2'] = 2 * N                     # Y2: chitobiose

    # Trimannosyl core — present in all mature N-glycans
    if n_N >= 2 and n_H >= 1:
        rem['+N2H'] = 2 * N + H               # +1 Man
    if n_N >= 2 and n_H >= 2:
        rem['+N2H2'] = 2 * N + 2 * H          # +2 Man
    if n_N >= 2 and n_H >= 3:
        rem['+N2H3'] = 2 * N + 3 * H          # trimannosyl core

    # =================================================================
    # High-mannose / paucimannose: linear mannose tree
    #   Man residues extend from the trimannosyl core on both arms.
    #   Sequential loss produces N2H3 → N2H4 → ... → N2H(n_hex).
    #   All intermediates are valid connected subgraphs.
    # =================================================================
    if gtype in ('high_mannose', 'paucimannose'):
        for h in range(4, n_H + 1):
            rem[f'+N2H{h}'] = 2 * N + h * H

    # =================================================================
    # Complex / hybrid: antenna-unit-based fragmentation
    #
    # Structure of each antenna on the trimannosyl core:
    #   Core-Man—GlcNAc—Gal(—NeuAc)
    #
    # In HCD, terminal residues are lost first (NeuAc, then Gal),
    # producing these valid intermediates above the core:
    #
    #   Core + n_ant × GlcNAc                (Gal all lost)
    #   Core + n_ant × GlcNAc + k × Gal     (k ≤ n_ant)
    #   Core + n_ant × GlcNAc + k × Gal + j × NeuAc  (j ≤ k)
    #
    # Biosynthetic rule: Gal count ≤ antenna GlcNAc count
    #   (Gal cannot exist without its preceding GlcNAc bridge)
    # =================================================================
    if gtype in ('complex', 'hybrid'):
        # Number of antennae = HexNAc beyond the chitobiose core
        n_ant = n_N - 2
        # Hex beyond the trimannosyl core = galactose on antennae
        n_gal = max(0, n_H - 3)

        # For hybrid: one arm is high-mannose (extra Man on core).
        # Treat extra Hex beyond (3 + n_ant) as mannose arm extension.
        extra_man = max(0, n_H - 3 - n_ant)

        # Hybrid mannose arm extension (added to core before antennae)
        if gtype == 'hybrid' and extra_man > 0:
            for m in range(1, extra_man + 1):
                h_count = 3 + m
                label = f'+N2H{h_count}'
                rem[label] = 2 * N + h_count * H

        # Antenna buildup: add GlcNAc one at a time, then Gal, then NeuAc
        for ant in range(1, n_ant + 1):
            core_n = 2 + ant     # chitobiose + antenna GlcNAc(s)
            core_h = 3           # trimannosyl core Man
            if gtype == 'hybrid':
                core_h += extra_man  # include hybrid mannose arm

            # Antenna GlcNAc only (all Gal stripped)
            label_base = f'+N{core_n}H{core_h}'
            rem[label_base] = core_n * N + core_h * H

            # Add Gal: up to min(ant, n_gal) — can't have more Gal than GlcNAc
            max_gal = min(ant, n_gal)
            for gal_k in range(1, max_gal + 1):
                h_total = core_h + gal_k
                label_g = f'+N{core_n}H{h_total}'
                rem[label_g] = core_n * N + h_total * H

                # Add NeuAc: up to min(gal_k, n_A) — NeuAc caps Gal
                for sia_j in range(1, min(gal_k, n_A) + 1):
                    label_s = f'+N{core_n}H{h_total}A{sia_j}' if sia_j > 1 else f'+N{core_n}H{h_total}A'
                    rem[label_s] = core_n * N + h_total * H + sia_j * A

    # =================================================================
    # Core fucose — only on the reducing-end GlcNAc (FUT8 product)
    # Add Fuc variant to chitobiose+ structures, not to monosaccharide
    # =================================================================
    if n_F >= 1:
        fuc_additions = {}
        for label, mass in list(rem.items()):
            # Core Fuc attaches to reducing-end GlcNAc → require N2+
            if label and mass >= 2 * N:
                fuc_additions[f'{label}F'] = mass + F
        rem.update(fuc_additions)

    # Intact glycopeptide
    rem['intact'] = total

    # Final filter: no partial glycan should exceed intact mass
    rem = {k: v for k, v in rem.items() if v <= total + 0.01}

    return rem


def generate_n_glycan_bions(glycan_comp):
    """
    Generate composition-specific glycan B-ions (oxonium fragments).

    Only returns B-ions whose monosaccharide content is a subset of the
    actual glycan composition.  Avoids false annotations (e.g. NeuAc
    oxonium when glycan has no sialic acid, or HexNAc-3Hex for complex
    glycans where a 3-Hex arm doesn't exist).

    Extended 2026-04-11 with the 17-branch diagnostic B-ion library from
    Shen et al. Nat Methods 2021 MOESM4 (StrucGP). Adds NeuGc, di-sialyl,
    di-fucosyl (Lewis-y), LacdiNAc, and sialyl-Lewis-x/a fragment ions
    when composition allows.

    Returns dict: {name: theoretical_mz} — all singly charged [M+H]+.
    """
    PROTON = 1.007276
    n_N = glycan_comp.get('HexNAc', 0)
    n_H = glycan_comp.get('Hex', 0)
    n_F = glycan_comp.get('Fuc', 0)
    n_A = glycan_comp.get('NeuAc', 0)
    n_G = glycan_comp.get('NeuGc', 0)

    gtype = _classify_glycan_type_simple(glycan_comp)

    # Monosaccharide residue masses
    _N = 203.0794    # HexNAc
    _H = 162.0528    # Hex
    _F = 146.0579    # Fuc (deoxyhexose)
    _A = 291.0954    # NeuAc
    _G = 307.0903    # NeuGc (NeuAc + 16 Da)

    ions = {}

    # --- Always present for N-glycans ---
    ions['HexNAc'] = _N + PROTON                    # 204.087  (B-1: terminal GlcNAc)
    ions['HexNAc-H2O'] = _N - 18.0106 + PROTON      # 186.076
    ions['HexNAc-2H2O'] = _N - 2*18.0106 + PROTON   # 168.066
    ions['HexNAc_138'] = 138.0550                     # cross-ring
    ions['Hex'] = _H + PROTON                         # 163.060
    ions['Hex-H2O'] = _H - 18.0106 + PROTON          # 145.050

    # --- Disaccharides (always: all N-glycans have HexNAc+Hex in core) ---
    ions['HexNAc-Hex'] = _N + _H + PROTON            # 366.140  (B-2: LacNAc)
    ions['HexNAc-Hex-H2O'] = _N + _H - 18.0106 + PROTON  # 348.129

    # --- Composition-dependent monosaccharides ---
    if n_F >= 1:
        ions['Fuc'] = _F + PROTON                     # 147.065
    if n_A >= 1:
        ions['NeuAc'] = _A + PROTON                   # 292.103  (B-3: NeuAc oxonium)
        ions['NeuAc-H2O'] = _A - 18.0106 + PROTON    # 274.092
    if n_G >= 1:
        ions['NeuGc'] = _G + PROTON                   # 308.098  (B-4: NeuGc oxonium)
        ions['NeuGc-H2O'] = _G - 18.0106 + PROTON    # 290.087

    # --- Composition-dependent disaccharides ---
    if n_A >= 1:
        ions['HexNAc-NeuAc'] = _N + _A + PROTON      # 495.182  (B-3 fragment)
    if n_G >= 1:
        ions['HexNAc-NeuGc'] = _N + _G + PROTON      # 511.177  (B-4 fragment)
        ions['Hex-NeuGc'] = _H + _G + PROTON         # 470.151  (B-10 fragment)
    if n_F >= 1:
        ions['HexNAc-Fuc'] = _N + _F + PROTON         # 350.145  (core Fuc or Lewis fragment)

    # --- Trisaccharides+ (glycan-type dependent) ---
    if gtype in ('high_mannose', 'paucimannose', 'hybrid'):
        # Mannose arm fragments: HexNAc + n×Hex
        if n_H >= 2:
            ions['HexNAc-2Hex'] = _N + 2*_H + PROTON   # 528.192
        if n_H >= 3:
            ions['HexNAc-3Hex'] = _N + 3*_H + PROTON   # 690.245
        if n_H >= 4:
            ions['HexNAc-4Hex'] = _N + 4*_H + PROTON   # 852.298

    if gtype in ('complex', 'hybrid'):
        # Sialyl-LacNAc fragment: NeuAc-Hex-HexNAc (B-3 full)
        if n_A >= 1 and n_H >= 4:  # Gal exists on antenna
            ions['HexNAc-Hex-NeuAc'] = _N + _H + _A + PROTON  # 657.235
        # NeuGc-LacNAc fragment (B-4 full)
        if n_G >= 1 and n_H >= 4:
            ions['HexNAc-Hex-NeuGc'] = _N + _H + _G + PROTON  # 673.230

    # Fucosylated trisaccharide (Lewis epitope): HexNAc-Hex-Fuc
    if n_F >= 1 and gtype in ('complex', 'hybrid'):
        ions['HexNAc-Hex-Fuc'] = _N + _H + _F + PROTON  # 512.197 (B-5: Lewis-x/a)

    # ========================================================================
    # Branch library extensions (Shen 2021 MOESM4)
    # Composition-conditional additions for complex/hybrid glycans
    # ========================================================================

    if gtype in ('complex', 'hybrid'):
        # B-6: Disialyl-LacNAc — NeuAc-NeuAc fragment (m/z 454.156)
        if n_A >= 2:
            ions['NeuAc-NeuAc'] = 2*_A - 18.0106 + PROTON   # 583.18 (NeuAc dimer with H2O loss)
            # Actual disialyl marker is NeuAc-Hex-HexNAc with NeuAc loss = 454.156
            # This is generated as a specific B-ion: 2*NeuAc - H2O - NeuAc + ...
            # Use the empirical value from MOESM4
            ions['Disialyl-marker_454'] = 454.156           # B-6 specific

        # B-8: Sialyl-Lewis-x/a — NeuAc-Fuc-Hex-HexNAc (m/z 803.293)
        if n_A >= 1 and n_F >= 1:
            ions['NeuAc-Fuc-Hex-HexNAc'] = _N + _H + _F + _A + PROTON  # 803.293

        # B-11: NeuGc-Lewis — NeuGc-Fuc-Hex-HexNAc (m/z 819.288)
        if n_G >= 1 and n_F >= 1:
            ions['NeuGc-Fuc-Hex-HexNAc'] = _N + _H + _F + _G + PROTON  # 819.288

        # B-12: Lewis-y di-fucosylated — HexNAc-Hex-2Fuc (m/z 658.272)
        if n_F >= 2:
            ions['HexNAc-Hex-2Fuc'] = _N + _H + 2*_F + PROTON  # 658.272

        # B-13: LacdiNAc — HexNAc-HexNAc dimer (m/z 407.166)
        # Only when extra HexNAc beyond core+antenna possible (n_N > 4 typically)
        if n_N >= 4:
            ions['HexNAc-HexNAc'] = 2*_N + PROTON           # 407.166 (B-13 LDN marker)
            # B-13 also: HexNAc-Hex-HexNAc (m/z 569.219)
            ions['HexNAc-Hex-HexNAc'] = 2*_N + _H + PROTON  # 569.219 (LDN or bisect)

        # B-14: LacdiNAc + Fuc (m/z 553.214)
        if n_N >= 4 and n_F >= 1:
            ions['HexNAc-HexNAc-Fuc'] = 2*_N + _F + PROTON  # 553.214

        # B-15: NeuAc-LacdiNAc (m/z 698.262)
        if n_N >= 4 and n_A >= 1:
            ions['NeuAc-HexNAc-HexNAc'] = 2*_N + _A + PROTON  # 698.262

        # B-7, B-9: Extended di-sialyl fragments
        if n_A >= 2:
            ions['NeuAc-NeuAc-Hex-HexNAc'] = _N + _H + 2*_A + PROTON  # 948.328
        if n_A >= 2 and n_F >= 1:
            ions['NeuAc-NeuAc-Fuc-Hex-HexNAc'] = _N + _H + _F + 2*_A + PROTON  # 1239.435
            ions['NeuAc-Fuc-Hex_intermediate'] = _H + _F + _A - 18.0106 + PROTON  # 583.218

    return ions
