# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Spectrum Annotator Ddzby is a Python package for MS/MS spectrum annotation of glycopeptides and crosslinked peptides. It calculates theoretical fragment ions, matches them against experimental spectra, and generates publication-ready IPSA-style plots.

## Common Commands

```bash
# Install for development
pip install -e .

# Install with Streamlit web UI
pip install -e ".[streamlit]"

# Install dev dependencies (pytest, black, flake8)
pip install -e ".[dev]"

# Run tests
pytest tests/
pytest tests/test_fragment_calculator.py -v  # single test file

# Format code (line-length=100)
black .

# Lint
flake8

# Run web interface
streamlit run app.py
```

## Architecture

The package follows a three-layer pipeline: **calculation → matching → visualization**.

```
spectrum_annotator_ddzby/
├── fragment_calculator.py  # Core: calculates theoretical m/z for all ion types
├── annotator.py           # Matches peaks, calculates FMR, generates plots
├── glycan_library.py      # Glycan/crosslinker databases and utilities
├── spectrum_reader.py     # mzML file parsing via pyteomics
└── __init__.py           # Public API exports
```

**Data Flow:**
1. `FragmentCalculator` takes peptide sequence + modifications → produces `TheoreticalIon` objects
2. `match_peaks()` compares theoretical ions to experimental m/z → produces `MatchedIon` objects
3. `SpectrumAnnotator` generates annotated PDF plots with statistics and FMR

**Key Data Classes:**
- `Modification`: position, residue, mass, name
- `TheoreticalIon`: ion type (b/y/c/z/Y/oxonium), number, charge, m/z, sequence
- `MatchedIon`: extends TheoreticalIon with experimental m/z, intensity, ppm error
- `GlycanComposition`: name, monosaccharide composition, mass
- `Crosslinker`: name, cleavable flag, stub masses (A/T/S for MS-cleavable)

## Key Technical Details

**Fragmentation Methods:** EThcD (default), HCD, CID, ETD - each enables different ion types

**Ion Types:**
- b/y: HCD backbone ions
- c/z: ETD backbone ions
- Y: glycopeptide Y ion series (Y0 = peptide only, up to Y(intact))
- Oxonium: glycan diagnostic ions (low m/z region)
- XL-stub: crosslinker stub ions (DSSO/DSBSO)

**Color Scheme (IPSA-compatible):**
- b: #0d75bc (blue), y: #be202d (red)
- c: #07a14a (green), z: #f79420 (orange)
- Y: #8491B4 (purple-gray), Oxonium: #7E6148 (brown)

**False Match Rate:** Uses spectrum shifting method (Schulte et al., Anal. Chem. 2025) - shifts m/z values to estimate random matching rate

**Glycan Libraries:**
- `O_GLYCAN_COMPOSITIONS`: 12+ human O-glycan structures (Core 1-4, sialylated variants)
- `N_GLYCAN_COMPOSITIONS`: High-mannose, complex, hybrid N-glycans
- `CROSSLINKERS`: DSSO, DSBSO (MS-cleavable), BS3, DSS (non-cleavable)

## Web App (app.py)

Streamlit interface with tabs: Annotate Spectrum, Glycan Library, Crosslinker Library, About. Settings sidebar controls tolerance, fragmentation method, and max charge state.
