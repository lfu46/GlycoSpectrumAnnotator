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


# Common O-glycan compositions
O_GLYCAN_COMPOSITIONS = {
    # Simple O-GlcNAc/O-GalNAc
    'O-GlcNAc': GlycanComposition('O-GlcNAc', {'HexNAc': 1}, 203.0794, 'O-GlcNAc'),
    'O-GalNAc': GlycanComposition('O-GalNAc', {'HexNAc': 1}, 203.0794, 'O-GalNAc'),

    # Core 1 (T antigen): GalNAc-Gal
    'Core1': GlycanComposition('Core1', {'HexNAc': 1, 'Hex': 1}, 365.1322, 'O-glycan'),

    # Core 2: GalNAc(-Gal)(-GlcNAc)
    'Core2': GlycanComposition('Core2', {'HexNAc': 2, 'Hex': 1}, 568.2116, 'O-glycan'),

    # Sialylated Core 1
    'SialylCore1': GlycanComposition('Sialyl-Core1', {'HexNAc': 1, 'Hex': 1, 'NeuAc': 1}, 656.2276, 'O-glycan'),
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
    'HexNAc-Hex-Hex': 528.1923,      # HexNAc + 2Hex

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
    'Hex': 163.0601,
    'Hex-H2O': 145.0495,
    'NeuAc': 292.1027,
    'NeuAc-H2O': 274.0921,
    'Fuc': 147.0652,
    'HexNAc-Hex': 366.1395,
    'HexNAc-Hex-H2O': 348.1289,
    'HexNAc-NeuAc': 495.1821,
    'HexNAc-Hex-NeuAc': 657.2349,
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
