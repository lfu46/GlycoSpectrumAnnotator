"""
Spectrum Annotator Ddzby - Universal MS/MS Spectrum Annotation

A Python package for annotating mass spectrometry spectra of glycopeptides
and crosslinked peptides, with support for multiple fragmentation methods.

Supports:
- O-glycans (Core 1-4, sialylated, phosphorylated, sulfated)
- N-glycans (high-mannose, complex, hybrid)
- Crosslinked peptides (DSSO, DSBSO, BS3)
- EThcD, HCD, CID fragmentation
- Y ion series and oxonium ions
- False match rate calculation
"""

__version__ = "0.3.0"
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

from .ms1_window import (
    plot_ms1_isolation_window,
    MS1WindowResult,
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
    # N-glycan HCD fragment generators (composition-aware remainder b/y + B-ions)
    generate_n_glycan_remainders,
    generate_n_glycan_bions,
    # Utilities
    identify_glycan_from_mass,
    get_glycan_mass,
    parse_proforma_glycan,
    # Crosslinker support
    Crosslinker,
    CROSSLINKERS,
    DSSO,
    DSBSO,
    BS3,
    DSS,
    generate_crosslink_fragments,
    identify_crosslink_stubs,
    calculate_crosslink_fmr,
    parse_proforma_crosslink,
)

__all__ = [
    # Classes
    "FragmentCalculator",
    "SpectrumAnnotator",
    "TheoreticalIon",
    "MatchedIon",
    "FalseMatchRate",
    "GlycanComposition",
    "Crosslinker",
    # Functions
    "match_peaks",
    "calculate_false_match_rate",
    "calculate_annotation_statistics",
    "parse_modifications_from_string",
    "annotate_spectra_batch",
    "plot_ms1_isolation_window",
    "MS1WindowResult",
    "generate_y_ion_series",
    "generate_n_glycan_y_ions",
    "generate_n_glycan_remainders",
    "generate_n_glycan_bions",
    "identify_glycan_from_mass",
    "get_glycan_mass",
    "parse_proforma_glycan",
    # Crosslinker functions
    "generate_crosslink_fragments",
    "identify_crosslink_stubs",
    "calculate_crosslink_fmr",
    "parse_proforma_crosslink",
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
    # Constants - Crosslinkers
    "CROSSLINKERS",
    "DSSO",
    "DSBSO",
    "BS3",
    "DSS",
    # Visualization
    "ION_COLORS",
]
