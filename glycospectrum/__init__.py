"""
GlycoSpectrumAnnotator - MS/MS Spectrum Annotation for Glycopeptides

A Python package for annotating mass spectrometry spectra of glycopeptides,
with support for EThcD fragmentation, Y ion series, oxonium ions, and
false match rate calculation.

Supports:
- O-GlcNAc / O-GalNAc (simple O-glycans)
- N-glycans (high-mannose, complex, hybrid)
- Extended Y ion series for complex glycans
- Comprehensive oxonium ion library
"""

__version__ = "0.1.0"
__author__ = "Longping Fu"

from .fragment_calculator import (
    FragmentCalculator,
    TheoreticalIon,
    MatchedIon,
    FalseMatchRate,
    match_peaks,
    calculate_false_match_rate,
    calculate_annotation_statistics,
    parse_modifications_from_string,
    OXONIUM_IONS,
    OXONIUM_IONS_ALL,
    AA_MASSES,
    MOD_MASSES,
)

from .annotator import (
    SpectrumAnnotator,
    annotate_spectra_batch,
    ION_COLORS,
)

from .glycan_library import (
    # Monosaccharide masses
    MONOSACCHARIDE_MASSES,
    # Oxonium ion sets
    OXONIUM_IONS_EXTENDED,
    OXONIUM_IONS_O_GLCNAC,
    OXONIUM_IONS_N_GLYCAN,
    # Glycan compositions
    GlycanComposition,
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    # Y ion generation
    generate_y_ion_series,
    generate_n_glycan_y_ions,
    # Utilities
    identify_glycan_from_mass,
    get_glycan_mass,
    parse_proforma_glycan,
)

__all__ = [
    # Classes
    "FragmentCalculator",
    "SpectrumAnnotator",
    "TheoreticalIon",
    "MatchedIon",
    "FalseMatchRate",
    "GlycanComposition",
    # Functions
    "match_peaks",
    "calculate_false_match_rate",
    "calculate_annotation_statistics",
    "parse_modifications_from_string",
    "annotate_spectra_batch",
    "generate_y_ion_series",
    "generate_n_glycan_y_ions",
    "identify_glycan_from_mass",
    "get_glycan_mass",
    "parse_proforma_glycan",
    # Constants - Masses
    "MONOSACCHARIDE_MASSES",
    "AA_MASSES",
    "MOD_MASSES",
    # Constants - Oxonium ions
    "OXONIUM_IONS",
    "OXONIUM_IONS_ALL",
    "OXONIUM_IONS_EXTENDED",
    "OXONIUM_IONS_O_GLCNAC",
    "OXONIUM_IONS_N_GLYCAN",
    # Constants - Glycan compositions
    "O_GLYCAN_COMPOSITIONS",
    "N_GLYCAN_COMPOSITIONS",
    # Visualization
    "ION_COLORS",
]
