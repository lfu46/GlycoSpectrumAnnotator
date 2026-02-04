"""
GlycoSpectrumAnnotator - MS/MS Spectrum Annotation for Glycopeptides

A Python package for annotating mass spectrometry spectra of glycopeptides,
with support for EThcD fragmentation, Y ion series, oxonium ions, and
false match rate calculation.
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
    AA_MASSES,
    MOD_MASSES,
)

from .annotator import (
    SpectrumAnnotator,
    annotate_spectra_batch,
    ION_COLORS,
)

__all__ = [
    # Classes
    "FragmentCalculator",
    "SpectrumAnnotator",
    "TheoreticalIon",
    "MatchedIon",
    "FalseMatchRate",
    # Functions
    "match_peaks",
    "calculate_false_match_rate",
    "calculate_annotation_statistics",
    "parse_modifications_from_string",
    "annotate_spectra_batch",
    # Constants
    "OXONIUM_IONS",
    "AA_MASSES",
    "MOD_MASSES",
    "ION_COLORS",
]
