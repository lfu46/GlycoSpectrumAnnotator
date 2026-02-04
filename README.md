# GlycoSpectrumAnnotator

A Python tool for annotating MS/MS spectra of glycopeptides, optimized for EThcD fragmentation with false match rate (FMR) calculation.

## Features

- **Glycopeptide-focused**: Optimized for O-GlcNAc and O-GalNAc modifications
- **EThcD support**: Full c/z ion series for glycan localization
- **Y ion annotation**: Intact glycopeptide and glycan loss ions (Y0/Y1)
- **Oxonium ions**: Diagnostic glycan fragment ions
- **False Match Rate**: Quality assessment using spectrum shifting method (Schulte et al., Anal. Chem. 2025)
- **Publication-ready**: IPSA-style colored PDF output
- **Batch processing**: Annotate multiple spectra with statistics export

## Installation

```bash
# Clone the repository
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator

# Install in development mode
pip install -e .

# Or install dependencies only
pip install -r requirements.txt
```

## Quick Start

### Single Spectrum Annotation

```python
from glycospectrum import SpectrumAnnotator, parse_modifications_from_string
import pandas as pd

# Load your spectrum data
spectrum = pd.read_csv("spectrum.csv")
exp_mz = spectrum['mz'].values
exp_intensity = spectrum['intensity'].values

# Define peptide and modifications
peptide = "AGYSQGATQYTQAQQTR"
modifications = parse_modifications_from_string("N-term(229.1629),4S(528.2859)")

# Create annotator
annotator = SpectrumAnnotator(
    peptide=peptide,
    modifications=modifications,
    precursor_charge=3,
    precursor_mz=789.1234,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,
    gene="GENE1",
    site_index="GENE1_S4"
)

# Generate annotated spectrum
fig = annotator.plot(output_path="annotated_spectrum.pdf")
```

### Batch Annotation

```python
from glycospectrum import annotate_spectra_batch

output_files, stats_df = annotate_spectra_batch(
    summary_file="spectra_summary.csv",
    spectra_dir="spectra/",
    output_dir="annotated/",
    tolerance_ppm=20.0,
    save_statistics=True
)

# View statistics
print(stats_df[['site_index', 'sequence_coverage', 'fmr_peaks']])
```

## Ion Types and Colors

| Ion Type | Color | Description |
|----------|-------|-------------|
| b | Blue | N-terminal HCD ions |
| y | Red | C-terminal HCD ions |
| c | Green | N-terminal ETD ions |
| z | Orange | C-terminal ETD ions |
| Y0 | Purple | Peptide with glycan loss |
| Y1 | Gray-purple | Intact glycopeptide |
| Oxonium | Brown | Glycan diagnostic ions |

## Output Statistics

The annotator provides comprehensive statistics:

- **Sequence Coverage**: Fraction of peptide bonds with fragment ions
- **FMR (peaks)**: False match rate based on peak count
- **FMR (intensity)**: False match rate based on intensity
- **Peaks Annotated**: Fraction of experimental peaks matched
- **Intensity Annotated**: Fraction of total intensity matched

## Supported Modifications

Default modifications include:
- TMT6plex (229.1629 Da)
- HexNAc / O-GlcNAc (203.0794 Da)
- HexNAc + TMT (528.2859 Da)
- Carbamidomethyl (57.02146 Da)
- Oxidation (15.9949 Da)

## File Formats

**Input**: CSV files with columns `mz` and `intensity`

**Output**:
- PDF annotated spectra (publication-ready)
- CSV annotation statistics

## References

- False Match Rate calculation: Schulte et al., Analytical Chemistry (2025)
- IPSA color scheme: [Interactive Peptide Spectral Annotator](https://github.com/coongroup/IPSA)

## Citation

If you use this tool in your research, please cite:

```
Fu, L. (2025). GlycoSpectrumAnnotator: MS/MS spectrum annotation for glycopeptides.
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
