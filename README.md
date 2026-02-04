# Spectrum Annotator Ddzby

Universal MS/MS spectrum annotation tool for glycopeptides and crosslinked peptides.

[![PyPI version](https://badge.fury.io/py/spectrum-annotator-ddzby.svg)](https://badge.fury.io/py/spectrum-annotator-ddzby)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

### Glycopeptide Annotation
- **O-glycans**: Core 1-4 structures, sialylated, phosphorylated, sulfated variants
- **N-glycans**: High-mannose, complex, hybrid structures
- **Y ion series**: Complete glycan ladder (Y0, Y1, Y2, ..., Y(core), Y(intact))
- **Oxonium ions**: Comprehensive diagnostic ion library

### Crosslinked Peptide Annotation
- **MS-cleavable**: DSSO, DSBSO with stub mass identification
- **Non-cleavable**: BS3, DSS
- **ProForma notation**: Full support for crosslink annotation

### Quality Assessment
- **False Match Rate**: Spectrum shifting method (Schulte et al., Anal. Chem. 2025)
- **Annotation statistics**: Coverage, matched peaks, intensity

## Installation

### From PyPI (Recommended)

```bash
pip install spectrum-annotator-ddzby
```

### From Source

```bash
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
pip install -e .
```

### With Streamlit Web App

```bash
pip install spectrum-annotator-ddzby[streamlit]
```

## Quick Start

### Python API

```python
from spectrum_annotator_ddzby import (
    FragmentCalculator,
    SpectrumAnnotator,
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    CROSSLINKERS,
    generate_crosslink_fragments,
)

# Glycopeptide annotation
calc = FragmentCalculator(
    peptide="PEPTIDEK",
    modifications={"5": "HexNAc"},
    precursor_charge=3,
    glycan_type="O-GlcNAc"
)
fragments = calc.calculate_all_fragments()

# Access glycan library
print(O_GLYCAN_COMPOSITIONS['Sialyl-T'].mass)  # 656.2276 Da

# Crosslinked peptide fragments
xl_fragments = generate_crosslink_fragments(
    peptide1_mass=1000.0,
    peptide2_mass=1200.0,
    crosslinker=CROSSLINKERS['DSBSO'],
    precursor_charge=4
)
```

### Web Interface (Streamlit)

```bash
# Run the web app
streamlit run app.py
```

Or visit the hosted version at: [Coming Soon]

## Supported Glycan Compositions

### O-Glycans (12 common human structures)
| Name | Composition | Mass (Da) |
|------|-------------|-----------|
| Tn (O-GlcNAc) | HexNAc(1) | 203.08 |
| T-antigen | HexNAc(1)Hex(1) | 365.13 |
| Sialyl-Tn | HexNAc(1)NeuAc(1) | 494.17 |
| Sialyl-T | HexNAc(1)Hex(1)NeuAc(1) | 656.23 |
| Core 2 | HexNAc(2)Hex(1) | 568.21 |
| ... | ... | ... |

### N-Glycans
| Name | Composition | Mass (Da) |
|------|-------------|-----------|
| Man5 | HexNAc(2)Hex(5) | 1216.42 |
| Man9 | HexNAc(2)Hex(9) | 1864.63 |
| A2G2F | HexNAc(4)Hex(5)Fuc(1) | 1768.64 |
| A2G2FS2 | HexNAc(4)Hex(5)Fuc(1)NeuAc(2) | 2350.83 |

## Supported Crosslinkers

| Name | Type | Spacer | Stub Masses |
|------|------|--------|-------------|
| DSSO | MS-cleavable | 10.1 Å | A: 54.01, T: 85.98, S: 103.99 |
| DSBSO | MS-cleavable | 12.5 Å | A: 54.01, T: 85.98, S: 103.99 |
| BS3 | Non-cleavable | 11.4 Å | - |
| DSS | Non-cleavable | 11.4 Å | - |

## Ion Types and Colors

| Ion Type | Color | Description |
|----------|-------|-------------|
| b | Blue | N-terminal HCD ions |
| y | Red | C-terminal HCD ions |
| c | Green | N-terminal ETD ions |
| z | Orange | C-terminal ETD ions |
| Y | Purple | Glycopeptide Y ions |
| Oxonium | Brown | Glycan diagnostic ions |
| XL-stub | Cyan | Crosslinker stubs |

## Output Statistics

- **Sequence Coverage**: Fraction of peptide bonds with fragment ions
- **FMR (peaks)**: False match rate based on peak count
- **FMR (intensity)**: False match rate based on intensity
- **Peaks Annotated**: Fraction of experimental peaks matched
- **Intensity Annotated**: Fraction of total intensity matched

## References

- Schulte, D. et al. "A Universal Spectrum Annotator for Complex Peptidoforms in Mass Spectrometry-Based Proteomics." *Analytical Chemistry* 2025, 97, 23120-23130.
- MSFragger Human O-glycan database

## Citation

If you use this tool in your research, please cite:

```
Fu, L. (2026). Spectrum Annotator Ddzby: Universal MS/MS spectrum annotation tool.
GitHub: https://github.com/lfu46/GlycoSpectrumAnnotator
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Author**: Longping Fu
- **Email**: lpfu46@gatech.edu
- **GitHub**: [@lfu46](https://github.com/lfu46)
