# Usage Guide

## Installation

```bash
pip install -e .
```

## Basic Usage

### Single Spectrum Annotation

```python
from glycospectrum import SpectrumAnnotator, parse_modifications_from_string
import pandas as pd

# Load spectrum
spectrum = pd.read_csv("spectrum.csv")

# Create annotator
annotator = SpectrumAnnotator(
    peptide="PEPTIDEK",
    modifications=parse_modifications_from_string("N-term(229.1629),5S(528.2859)"),
    precursor_charge=3,
    precursor_mz=500.0,
    exp_mz=spectrum['mz'].values,
    exp_intensity=spectrum['intensity'].values,
)

# Generate plot
fig = annotator.plot(output_path="annotated.pdf")
```

### Batch Processing

```python
from glycospectrum import annotate_spectra_batch

files, stats = annotate_spectra_batch(
    summary_file="summary.csv",
    spectra_dir="spectra/",
    output_dir="output/",
)
```

## Configuration

See `config/default_config.yaml` for all available options.
