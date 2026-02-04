# User Manual: Spectrum Annotator Ddzby

A comprehensive guide for installing and using Spectrum Annotator Ddzby on Windows, macOS, and Linux.

---

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation](#installation)
   - [Windows](#windows-installation)
   - [macOS](#macos-installation)
   - [Linux](#linux-installation)
3. [Quick Start](#quick-start)
4. [Python API Usage](#python-api-usage)
   - [Calculating Theoretical Fragments](#calculating-theoretical-fragments)
   - [Matching Peaks](#matching-peaks)
   - [Creating Annotated Spectrum Plots](#creating-annotated-spectrum-plots)
   - [Batch Processing](#batch-processing)
   - [Working with Glycans](#working-with-glycans)
   - [Working with Crosslinkers](#working-with-crosslinkers)
5. [Web Interface (Streamlit)](#web-interface-streamlit)
6. [Input Formats](#input-formats)
7. [Output Formats](#output-formats)
8. [Configuration Options](#configuration-options)
9. [Troubleshooting](#troubleshooting)

---

## System Requirements

- **Python**: 3.9 or higher (tested on 3.9, 3.10, 3.11, 3.12)
- **Operating System**: Windows 10/11, macOS 10.15+, Linux (Ubuntu 20.04+, CentOS 7+)
- **Memory**: 4 GB RAM minimum, 8 GB recommended for batch processing
- **Disk Space**: 100 MB for installation

---

## Installation

### Windows Installation

#### Option 1: Using pip (Recommended)

1. **Install Python** (if not already installed):
   - Download Python from https://www.python.org/downloads/windows/
   - During installation, check "Add Python to PATH"
   - Verify installation by opening Command Prompt and running:
     ```cmd
     python --version
     ```

2. **Install the package**:
   ```cmd
   pip install spectrum-annotator-ddzby
   ```

3. **Install with web interface** (optional):
   ```cmd
   pip install spectrum-annotator-ddzby[streamlit]
   ```

#### Option 2: From Source

1. **Install Git** from https://git-scm.com/download/win

2. **Clone and install**:
   ```cmd
   git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
   cd GlycoSpectrumAnnotator
   pip install -e .
   ```

### macOS Installation

#### Option 1: Using pip (Recommended)

1. **Install Python** (if not using system Python):
   ```bash
   # Using Homebrew (recommended)
   brew install python@3.11

   # Or download from https://www.python.org/downloads/macos/
   ```

2. **Install the package**:
   ```bash
   pip3 install spectrum-annotator-ddzby
   ```

3. **Install with web interface** (optional):
   ```bash
   pip3 install "spectrum-annotator-ddzby[streamlit]"
   ```

#### Option 2: From Source

```bash
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
pip3 install -e .
```

### Linux Installation

#### Ubuntu/Debian

1. **Install Python and pip**:
   ```bash
   sudo apt update
   sudo apt install python3 python3-pip python3-venv
   ```

2. **Create virtual environment** (recommended):
   ```bash
   python3 -m venv glyco_env
   source glyco_env/bin/activate
   ```

3. **Install the package**:
   ```bash
   pip install spectrum-annotator-ddzby
   ```

#### CentOS/RHEL/Fedora

1. **Install Python**:
   ```bash
   # CentOS/RHEL 8+
   sudo dnf install python39 python39-pip

   # Fedora
   sudo dnf install python3 python3-pip
   ```

2. **Install the package**:
   ```bash
   pip3 install spectrum-annotator-ddzby
   ```

#### Conda Environment (All Platforms)

```bash
# Create and activate environment
conda create -n glyco python=3.11
conda activate glyco

# Install package
pip install spectrum-annotator-ddzby
```

---

## Quick Start

### Verify Installation

```python
from spectrum_annotator_ddzby import __version__
print(f"Spectrum Annotator Ddzby version: {__version__}")
```

### Basic Example

```python
from spectrum_annotator_ddzby import (
    FragmentCalculator,
    O_GLYCAN_COMPOSITIONS,
)

# Define a glycopeptide
peptide = "PEPTIDEK"
modifications = [
    {'position': 0, 'residue': 'N-term', 'mass': 229.1629},  # TMT
    {'position': 5, 'residue': 'T', 'mass': 203.0794},       # O-GlcNAc
]

# Calculate theoretical fragments
calc = FragmentCalculator(
    peptide=peptide,
    modifications=modifications,
    precursor_charge=3,
    max_fragment_charge=2
)

print(f"Precursor mass: {calc.precursor_mass:.4f} Da")
print(f"Precursor m/z: {calc.precursor_mz:.4f}")

# Get all theoretical ions
all_ions = calc.calculate_all_ions()
for ion_type, ions in all_ions.items():
    if ions:
        print(f"{ion_type}: {len(ions)} ions")
```

---

## Python API Usage

### Calculating Theoretical Fragments

The `FragmentCalculator` class generates theoretical m/z values for all ion types.

```python
from spectrum_annotator_ddzby import FragmentCalculator

# Initialize calculator
calc = FragmentCalculator(
    peptide="AGYSQGATQYTQAQQTR",
    modifications=[
        {'position': 0, 'residue': 'N-term', 'mass': 229.1629},   # TMT
        {'position': 4, 'residue': 'S', 'mass': 528.2859},        # HexNAc+TMT
    ],
    precursor_charge=3,
    max_fragment_charge=2,
    glycan_type='O-GlcNAc',          # 'auto', 'O-GlcNAc', 'O-GalNAc', 'N-glycan'
    use_extended_oxonium=False        # True for N-glycans
)

# Calculate specific ion types
b_ions = calc.calculate_b_ions()           # N-terminal HCD
y_ions = calc.calculate_y_ions()           # C-terminal HCD
c_ions = calc.calculate_c_ions()           # N-terminal ETD
z_ions = calc.calculate_z_ions()           # C-terminal ETD
Y_ions = calc.calculate_Y_ions()           # Glycopeptide Y ions
oxonium = calc.calculate_oxonium_ions()    # Glycan diagnostic ions

# Calculate all ions at once
all_ions = calc.calculate_all_ions(
    include_neutral_losses=True,
    neutral_loss_types=['H2O', 'NH3', 'HexNAc_TMT', 'HexNAc']
)

# Get flat list of all ions
all_ions_flat = calc.get_all_ions_flat()
print(f"Total theoretical ions: {len(all_ions_flat)}")
```

### Matching Peaks

Match theoretical ions against experimental spectrum data.

```python
import numpy as np
from spectrum_annotator_ddzby import (
    FragmentCalculator,
    match_peaks,
    calculate_false_match_rate,
    calculate_annotation_statistics,
)

# Example experimental data
exp_mz = np.array([204.087, 366.140, 500.25, 650.30, 800.40])
exp_intensity = np.array([1000, 800, 500, 300, 200])

# Calculate theoretical ions
calc = FragmentCalculator(
    peptide="PEPTIDEK",
    modifications=[{'position': 5, 'residue': 'T', 'mass': 203.0794}],
    precursor_charge=3
)
theoretical_ions = calc.get_all_ions_flat()

# Match peaks
matched_ions = match_peaks(
    theoretical_ions=theoretical_ions,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,           # Mass tolerance in ppm
    match_isotopes=True,          # Also match M+1, M+2 isotopes
    max_isotope=2                 # Maximum isotope offset
)

print(f"Matched ions: {len(matched_ions)}")
for ion in matched_ions[:5]:
    print(f"  {ion.annotation}: theo={ion.mz:.4f}, exp={ion.exp_mz:.4f}, "
          f"error={ion.mass_error_ppm:.1f} ppm")

# Calculate False Match Rate (spectrum shifting method)
fmr = calculate_false_match_rate(
    theoretical_ions=theoretical_ions,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,
    shift_range=25.0,             # Shift range in Th
    shift_step=1.0                # Shift step size
)

print(f"\nFalse Match Rate:")
print(f"  FMR (peaks): {fmr.fmr_peaks*100:.1f}%")
print(f"  FMR (intensity): {fmr.fmr_intensity*100:.1f}%")
print(f"  Matched peaks: {fmr.matched_peaks}")
print(f"  Avg random matches: {fmr.avg_random_peaks:.2f}")

# Calculate annotation statistics
stats = calculate_annotation_statistics(
    matched_ions=matched_ions,
    theoretical_ions=theoretical_ions,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    peptide_length=len("PEPTIDEK")
)

print(f"\nAnnotation Statistics:")
print(f"  Sequence coverage: {stats['sequence_coverage']*100:.1f}%")
print(f"  Peaks annotated: {stats['peaks_annotated']*100:.1f}%")
print(f"  Intensity annotated: {stats['intensity_annotated']*100:.1f}%")
```

### Creating Annotated Spectrum Plots

Generate publication-quality IPSA-style annotated spectra.

```python
import numpy as np
import pandas as pd
from spectrum_annotator_ddzby import SpectrumAnnotator

# Load spectrum data (CSV with 'mz' and 'intensity' columns)
# spectrum_df = pd.read_csv("spectrum.csv")
# exp_mz = spectrum_df['mz'].values
# exp_intensity = spectrum_df['intensity'].values

# Example data
exp_mz = np.array([204.087, 366.140, 500.25, 650.30, 800.40, 950.50])
exp_intensity = np.array([1000, 800, 500, 300, 200, 150])

# Create annotator
annotator = SpectrumAnnotator(
    peptide="AGYSQGATQYTQAQQTR",
    modifications=[
        {'position': 0, 'residue': 'N-term', 'mass': 229.1629},
        {'position': 4, 'residue': 'S', 'mass': 528.2859},
    ],
    precursor_charge=3,
    precursor_mz=800.5,                # Experimental precursor m/z
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,
    site_index="Q96KR1_S195",          # Site identifier
    gene="GENE1"                        # Gene name
)

# Generate plot
fig = annotator.plot(
    figsize=(7, 5),
    output_path="annotated_spectrum.pdf",   # Save to PDF
    show_error_plot=True,                    # Include mass error panel
    intensity_threshold_pct=1.0,             # Min intensity to label (% base peak)
    max_labels=50                            # Maximum labels on plot
)

# Access computed statistics
print(f"Sequence coverage: {annotator.annotation_stats['sequence_coverage']*100:.1f}%")
print(f"FMR (peaks): {annotator.false_match_rate.fmr_peaks*100:.1f}%")
```

### Batch Processing

Annotate multiple spectra from a summary file.

```python
from spectrum_annotator_ddzby import annotate_spectra_batch

# Batch annotate spectra
output_files, stats_df = annotate_spectra_batch(
    summary_file="extracted_spectra_summary.csv",    # Summary CSV file
    spectra_dir="spectra/",                          # Directory with spectrum CSVs
    output_dir="annotated_spectra/",                 # Output directory for PDFs
    n_spectra=100,                                   # Number to process (None=all)
    tolerance_ppm=20.0,
    save_statistics=True                             # Save statistics CSV
)

print(f"Generated {len(output_files)} annotated spectra")

# View statistics summary
if stats_df is not None:
    print(f"\nMean sequence coverage: {stats_df['sequence_coverage'].mean()*100:.1f}%")
    print(f"Mean FMR (peaks): {stats_df['fmr_peaks'].mean()*100:.1f}%")
```

### Working with Glycans

Access the glycan composition library and calculate glycan masses.

```python
from spectrum_annotator_ddzby import (
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    MONOSACCHARIDE_MASSES,
    GlycanComposition,
    get_glycan_mass,
    identify_glycan_from_mass,
    generate_y_ion_series,
    parse_proforma_glycan,
)

# Access O-glycan library
print("O-Glycan Library:")
for name, glycan in O_GLYCAN_COMPOSITIONS.items():
    print(f"  {name}: {glycan.mass:.4f} Da - {glycan.composition}")

# Access N-glycan library
print("\nN-Glycan Library:")
for name, glycan in list(N_GLYCAN_COMPOSITIONS.items())[:5]:
    print(f"  {name}: {glycan.mass:.4f} Da")

# Access monosaccharide masses
print(f"\nMonosaccharide Masses:")
print(f"  HexNAc: {MONOSACCHARIDE_MASSES['HexNAc']:.4f} Da")
print(f"  Hex: {MONOSACCHARIDE_MASSES['Hex']:.4f} Da")
print(f"  NeuAc: {MONOSACCHARIDE_MASSES['NeuAc']:.4f} Da")
print(f"  Fuc: {MONOSACCHARIDE_MASSES['Fuc']:.4f} Da")

# Create custom glycan composition
custom_glycan = GlycanComposition.from_string(
    "HexNAc2Hex5Fuc1",
    name="Custom N-glycan",
    glycan_type="N-glycan"
)
print(f"\nCustom glycan: {custom_glycan.mass:.4f} Da")

# Calculate glycan mass from composition
mass = get_glycan_mass({'HexNAc': 2, 'Hex': 5, 'Fuc': 1, 'NeuAc': 1})
print(f"Calculated mass: {mass:.4f} Da")

# Identify glycan from mass
identified = identify_glycan_from_mass(656.2276, tolerance=0.01)
print(f"\nIdentified glycan: {identified}")

# Generate Y ion series for a glycan
peptide_mass = 1500.0  # Peptide mass without glycan
glycan = O_GLYCAN_COMPOSITIONS['Sialyl-T']
y_ions = generate_y_ion_series(peptide_mass, glycan.composition, charges=[2, 3])
for ion in y_ions[:5]:
    print(f"  {ion['name']}: m/z={ion['mz']:.4f} (charge {ion['charge']}+)")
```

### Working with Crosslinkers

Calculate fragments for crosslinked peptides.

```python
from spectrum_annotator_ddzby import (
    CROSSLINKERS,
    DSSO,
    DSBSO,
    BS3,
    DSS,
    Crosslinker,
    generate_crosslink_fragments,
    identify_crosslink_stubs,
    parse_proforma_crosslink,
)

# Access crosslinker library
print("Crosslinker Library:")
for name, xl in CROSSLINKERS.items():
    print(f"  {name}: mass={xl.intact_mass:.4f} Da, cleavable={xl.cleavable}")
    if xl.cleavable and xl.stub_masses:
        for stub_name, stub_mass in xl.stub_masses.items():
            print(f"    Stub {stub_name}: {stub_mass:.4f} Da")

# Generate crosslink fragments for MS-cleavable crosslinker
xl_fragments = generate_crosslink_fragments(
    peptide1_mass=1000.0,
    peptide2_mass=1200.0,
    crosslinker=CROSSLINKERS['DSSO'],
    precursor_charge=4
)

print(f"\nCrosslink fragments for DSSO:")
for frag in xl_fragments[:10]:
    print(f"  {frag['annotation']}: m/z={frag['mz']:.4f}")

# Identify stub masses in spectrum
stub_matches = identify_crosslink_stubs(
    exp_mz=np.array([1054.01, 1085.98, 1103.99]),
    exp_intensity=np.array([1000, 800, 600]),
    peptide_mass=1000.0,
    crosslinker=DSSO,
    tolerance_ppm=20.0
)

for match in stub_matches:
    print(f"  Found stub: {match['stub']} at m/z={match['exp_mz']:.4f}")
```

---

## Web Interface (Streamlit)

### Starting the Web App

```bash
# Navigate to the project directory
cd GlycoSpectrumAnnotator

# Run the Streamlit app
streamlit run app.py
```

The app will open in your default browser at http://localhost:8501

### Web Interface Features

#### 1. Annotate Spectrum Tab

- **Input Section**:
  - Enter peptide sequence
  - Select analysis type (Glycopeptide, Crosslinked Peptide, or Linear Peptide)
  - For glycopeptides: select glycan type and composition
  - For crosslinked peptides: enter second peptide and select crosslinker
  - Set precursor charge and m/z
  - Paste peak list (m/z and intensity, space-separated)

- **Settings Sidebar**:
  - Mass tolerance (ppm or Da)
  - Fragmentation method (EThcD, HCD, CID, ETD)
  - Maximum fragment charge

- **Results**:
  - Interactive annotated spectrum plot
  - Matched ion highlighting
  - Annotation summary statistics

#### 2. Glycan Library Tab

- Browse O-glycan and N-glycan compositions
- View mass, composition, and type for each glycan
- **Glycan Mass Calculator**: Calculate mass from monosaccharide counts

#### 3. Crosslinker Library Tab

- Browse crosslinker definitions
- View intact mass, spacer length, and MS-cleavability
- Stub mass details for MS-cleavable crosslinkers

#### 4. About Tab

- Installation instructions
- Python API examples
- Citation information

---

## Input Formats

### Modification String Format (FragPipe style)

```
N-term(229.1629),4S(528.2859),19K(229.1629)
```

- `N-term(mass)`: N-terminal modification
- `C-term(mass)`: C-terminal modification
- `<position><residue>(mass)`: Residue modification (1-indexed)

### Peak List Format

Tab or space-separated m/z and intensity values, one peak per line:

```
204.087 1000
366.140 800
500.25  500
650.30  300
800.40  200
```

### Spectrum CSV Format

CSV file with `mz` and `intensity` columns:

```csv
mz,intensity
204.087,1000
366.140,800
500.25,500
650.30,300
800.40,200
```

### Batch Summary CSV Format

Required columns for batch processing:
- `Peptide`: Peptide sequence
- `Charge`: Precursor charge state
- `site_index`: Site identifier
- `Gene`: Gene name
- `spectrum_file`: Path to spectrum CSV
- `Assigned_Modifications` or `modifications_json`: Modification info
- `Calibrated_Observed_MZ` or `mzML_precursor_mz`: Precursor m/z

---

## Output Formats

### Annotated Spectrum PDF

Publication-ready PDF containing:
- Peptide sequence with fragmentation sites
- Sample and site information
- Sequence coverage and False Match Rate
- Annotated spectrum with color-coded ions
- Mass error scatter plot
- Legend with ion type colors

### Ion Type Color Scheme

| Ion Type | Color | Description |
|----------|-------|-------------|
| b | Blue (#0d75bc) | N-terminal HCD ions |
| y | Red (#be202d) | C-terminal HCD ions |
| c | Green (#07a14a) | N-terminal ETD ions |
| z | Orange (#f79420) | C-terminal ETD ions |
| Y0 | Purple (#9B59B6) | Peptide with glycan loss |
| Y1 | Purple-gray (#8491B4) | Intact glycopeptide |
| Oxonium | Brown (#7E6148) | Glycan diagnostic ions |
| Unassigned | Gray (#a6a6a6) | Unmatched peaks |

### Statistics CSV (Batch Processing)

Output columns:
- `site_index`, `gene`, `peptide`, `charge`, `scan_number`
- `sequence_coverage`, `sequence_coverage_bonds`
- `peaks_annotated`, `peaks_annotated_count`
- `intensity_annotated`
- `fragments_found`, `fragments_found_count`
- `fmr_peaks`, `fmr_intensity`
- `matched_peaks`, `matched_intensity`
- `avg_random_peaks`, `avg_random_intensity`

---

## Configuration Options

### FragmentCalculator Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `peptide` | str | required | Peptide sequence |
| `modifications` | list | required | List of modification dicts |
| `precursor_charge` | int | required | Precursor charge state |
| `max_fragment_charge` | int | 2 | Maximum fragment charge |
| `glycan_type` | str | 'auto' | Glycan type hint |
| `use_extended_oxonium` | bool | False | Use extended oxonium library |

### match_peaks Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tolerance_ppm` | float | 20.0 | Mass tolerance in ppm |
| `match_isotopes` | bool | True | Match isotope peaks |
| `max_isotope` | int | 2 | Maximum isotope offset |

### SpectrumAnnotator.plot Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `figsize` | tuple | (7, 5) | Figure size (width, height) |
| `output_path` | str | None | Path to save PDF |
| `show_error_plot` | bool | True | Include mass error panel |
| `intensity_threshold_pct` | float | 1.0 | Min intensity to label (%) |
| `max_labels` | int | 50 | Maximum labels on plot |

---

## Troubleshooting

### Installation Issues

**Problem**: `pip install` fails with permission error
```bash
# Solution: Use --user flag or virtual environment
pip install --user spectrum-annotator-ddzby

# Or create virtual environment
python -m venv env
source env/bin/activate  # Linux/macOS
env\Scripts\activate     # Windows
pip install spectrum-annotator-ddzby
```

**Problem**: Import error for pyteomics
```bash
# Solution: Install pyteomics separately
pip install pyteomics lxml
```

**Problem**: Matplotlib font warnings on Linux
```bash
# Solution: Install fonts
sudo apt install fonts-liberation  # Ubuntu/Debian
sudo dnf install liberation-fonts  # Fedora/RHEL
```

### Runtime Issues

**Problem**: Memory error during batch processing
- Solution: Process fewer spectra at a time using `n_spectra` parameter
- Solution: Close figures after saving with `plt.close(fig)`

**Problem**: No peaks matched
- Check mass tolerance (try increasing to 30-50 ppm)
- Verify modifications are correct
- Check precursor charge state
- Ensure m/z values are in the expected range

**Problem**: Streamlit app won't start
```bash
# Check if port 8501 is in use
lsof -i :8501  # Linux/macOS
netstat -an | findstr 8501  # Windows

# Use different port
streamlit run app.py --server.port 8502
```

### Contact and Support

- **Author**: Longping Fu
- **Email**: lpfu46@gatech.edu
- **GitHub Issues**: https://github.com/lfu46/GlycoSpectrumAnnotator/issues

---

## References

- Schulte, D. et al. "A Universal Spectrum Annotator for Complex Peptidoforms in Mass Spectrometry-Based Proteomics." *Analytical Chemistry* 2025, 97, 23120-23130.
- MSFragger Human O-glycan database
