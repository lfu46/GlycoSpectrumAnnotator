# User Manual: Spectrum Annotator Ddzby

A comprehensive guide for installing and using Spectrum Annotator Ddzby on Windows, macOS, and Linux.

---

## Table of Contents

1. [Overview: Two Ways to Use This Tool](#overview-two-ways-to-use-this-tool)
2. [System Requirements](#system-requirements)
3. [Installation](#installation)
   - [Windows](#windows-installation)
   - [macOS](#macos-installation)
   - [Linux](#linux-installation)
4. [Using the Web Interface (No Coding Required)](#using-the-web-interface-no-coding-required)
5. [Using the Python API (For Developers)](#using-the-python-api-for-developers)
6. [Input and Output Formats](#input-and-output-formats)
7. [Troubleshooting](#troubleshooting)
8. [Reference](#reference)

---

## Overview: Two Ways to Use This Tool

Spectrum Annotator Ddzby can be used in **two ways**:

### 1. Web Interface (Recommended for Most Users)
- **No programming required**
- Point-and-click interface in your web browser
- Interactive spectrum visualization
- Great for annotating individual spectra

### 2. Python API (For Developers and Batch Processing)
- Requires writing Python code
- Powerful for batch processing hundreds of spectra
- Can be integrated into analysis pipelines
- More configuration options

**Which should you choose?**
- If you just want to annotate a few spectra → Use the **Web Interface**
- If you need to process many spectra automatically → Use the **Python API**

---

## System Requirements

- **Python**: 3.9 or higher
- **Operating System**: Windows 10/11, macOS 10.15+, or Linux
- **Memory**: 4 GB RAM minimum
- **Disk Space**: 100 MB

---

## Installation

### Windows Installation

#### Step 1: Install Python

1. Download Python from https://www.python.org/downloads/windows/
2. Run the installer
3. **IMPORTANT**: Check the box "Add Python to PATH" at the bottom of the installer
4. Click "Install Now"

#### Step 2: Open Command Prompt

1. Press `Windows + R`
2. Type `cmd` and press Enter
3. A black terminal window will open

#### Step 3: Verify Python is Installed

Type this command and press Enter:
```cmd
python --version
```
You should see something like `Python 3.11.5`

#### Step 4: Install the Package

Type this command and press Enter:
```cmd
pip install "spectrum-annotator-ddzby[streamlit]"
```
Wait for the installation to complete (this may take a few minutes).

#### Step 5: Download the Web App

```cmd
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
```

Or download manually:
1. Go to https://github.com/lfu46/GlycoSpectrumAnnotator
2. Click the green "Code" button
3. Click "Download ZIP"
4. Extract the ZIP file to a folder (e.g., `C:\GlycoSpectrumAnnotator`)

#### Step 6: Start the Web Interface

Navigate to the folder and run:
```cmd
cd C:\GlycoSpectrumAnnotator
streamlit run app.py
```

Your web browser will automatically open to `http://localhost:8501`

---

### macOS Installation

#### Step 1: Open Terminal

1. Press `Cmd + Space` to open Spotlight
2. Type `Terminal` and press Enter

#### Step 2: Check if Python is Installed

```bash
python3 --version
```

If you see `command not found`, install Python:
```bash
# Option A: Using Homebrew (recommended)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install python@3.11

# Option B: Download from https://www.python.org/downloads/macos/
```

#### Step 3: Install the Package

```bash
pip3 install "spectrum-annotator-ddzby[streamlit]"
```

#### Step 4: Download the Web App

```bash
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
```

Or download manually from https://github.com/lfu46/GlycoSpectrumAnnotator

#### Step 5: Start the Web Interface

```bash
cd ~/GlycoSpectrumAnnotator
streamlit run app.py
```

Your web browser will automatically open to `http://localhost:8501`

---

### Linux Installation

#### Ubuntu/Debian

```bash
# Install Python and pip
sudo apt update
sudo apt install python3 python3-pip git

# Install the package
pip3 install "spectrum-annotator-ddzby[streamlit]"

# Download and run
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
streamlit run app.py
```

#### CentOS/RHEL/Fedora

```bash
# Install Python
sudo dnf install python3 python3-pip git

# Install the package
pip3 install "spectrum-annotator-ddzby[streamlit]"

# Download and run
git clone https://github.com/lfu46/GlycoSpectrumAnnotator.git
cd GlycoSpectrumAnnotator
streamlit run app.py
```

---

## Using the Web Interface (No Coding Required)

The web interface is the easiest way to use Spectrum Annotator Ddzby.

### Starting the Web App

**Every time you want to use the tool**, open your terminal and run:

**Windows (Command Prompt):**
```cmd
cd C:\GlycoSpectrumAnnotator
streamlit run app.py
```

**macOS/Linux (Terminal):**
```bash
cd ~/GlycoSpectrumAnnotator
streamlit run app.py
```

The web browser will open automatically. If it doesn't, go to: http://localhost:8501

### Step-by-Step: Annotating a Spectrum

#### 1. Choose Your Analysis Type

In the left sidebar, select one of:
- **Glycopeptide** - For O-GlcNAc, O-GalNAc, or N-glycan modified peptides
- **Crosslinked Peptide** - For DSSO, DSBSO, BS3, or DSS crosslinked peptides
- **Linear Peptide** - For unmodified or simple modifications

#### 2. Set Mass Tolerance

In the sidebar under "Mass Tolerance":
- **ppm** (recommended): Use 10-20 ppm for high-resolution data
- **Da**: Use 0.02-0.05 Da for unit resolution data

#### 3. Enter Your Peptide Sequence

In the main panel, enter your peptide sequence using single-letter amino acid codes:
```
PEPTIDEK
```

#### 4. Configure Modifications (for Glycopeptides)

- Select **O-glycan** or **N-glycan**
- Choose the glycan composition from the dropdown (e.g., "Sialyl-T", "Man5")
- Enter the **Glycosylation Site** (position number, 1-indexed)

#### 5. Enter Precursor Information

- **Charge State**: The precursor charge (e.g., 2, 3, 4)
- **Precursor m/z**: Optional, for display purposes

#### 6. Paste Your Spectrum Data

In the "Spectrum Data" text box, paste your peak list:
```
204.087 1000
366.140 800
500.25 500
650.30 300
```

Format: `m/z` `intensity` (space or tab separated, one peak per line)

You can export this from:
- **Thermo Xcalibur**: File → Export → Peak List
- **Waters MassLynx**: File → Export → Text
- **Bruker DataAnalysis**: File → Export → Peak List

#### 7. Click "Annotate Spectrum"

The annotated spectrum will appear showing:
- Matched peaks highlighted in color
- Ion annotations (b, y, c, z, Y, oxonium)
- Summary statistics

### Using the Glycan Library

Click the **"Glycan Library"** tab to:
- Browse all O-glycan and N-glycan compositions
- View mass and composition for each glycan
- Use the **Glycan Mass Calculator** to calculate custom glycan masses

### Using the Crosslinker Library

Click the **"Crosslinker Library"** tab to:
- View all supported crosslinkers (DSSO, DSBSO, BS3, DSS)
- See stub masses for MS-cleavable crosslinkers

### Stopping the Web App

To stop the web app:
1. Go back to your terminal
2. Press `Ctrl + C`

---

## Using the Python API (For Developers)

For batch processing or integration into pipelines, use the Python API.

### Quick Start: Verify Installation

Open a Python interpreter or create a new Python script:

**In Terminal:**
```bash
python3
```

**Then type:**
```python
from spectrum_annotator_ddzby import __version__
print(f"Version: {__version__}")
```

### Example 1: Calculate Theoretical Fragments

Create a file called `example1.py`:

```python
from spectrum_annotator_ddzby import FragmentCalculator

# Define peptide and modifications
peptide = "PEPTIDEK"
modifications = [
    {'position': 0, 'residue': 'N-term', 'mass': 229.1629},  # TMT
    {'position': 5, 'residue': 'T', 'mass': 203.0794},       # O-GlcNAc
]

# Create calculator
calc = FragmentCalculator(
    peptide=peptide,
    modifications=modifications,
    precursor_charge=3
)

# Print results
print(f"Peptide: {peptide}")
print(f"Precursor mass: {calc.precursor_mass:.4f} Da")
print(f"Precursor m/z: {calc.precursor_mz:.4f}")

# Get all theoretical ions
all_ions = calc.calculate_all_ions()
for ion_type, ions in all_ions.items():
    if ions:
        print(f"  {ion_type}: {len(ions)} ions")
```

**Run it:**
```bash
python3 example1.py
```

### Example 2: Annotate a Spectrum and Create PDF

Create a file called `example2.py`:

```python
import numpy as np
from spectrum_annotator_ddzby import SpectrumAnnotator

# Your experimental spectrum data
exp_mz = np.array([204.087, 366.140, 500.25, 650.30, 800.40, 950.50])
exp_intensity = np.array([1000, 800, 500, 300, 200, 150])

# Create annotator
annotator = SpectrumAnnotator(
    peptide="AGYSQGATQYTQAQQTR",
    modifications=[
        {'position': 0, 'residue': 'N-term', 'mass': 229.1629},  # TMT
        {'position': 4, 'residue': 'S', 'mass': 528.2859},       # HexNAc+TMT
    ],
    precursor_charge=3,
    precursor_mz=800.5,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,
    site_index="PROTEIN_S4",
    gene="GENE1"
)

# Generate and save PDF
fig = annotator.plot(output_path="my_annotated_spectrum.pdf")
print("Saved: my_annotated_spectrum.pdf")

# Print statistics
print(f"Sequence coverage: {annotator.annotation_stats['sequence_coverage']*100:.1f}%")
print(f"False Match Rate: {annotator.false_match_rate.fmr_peaks*100:.1f}%")
```

**Run it:**
```bash
python3 example2.py
```

### Example 3: Load Spectrum from CSV File

If your spectrum is saved as a CSV file with columns `mz` and `intensity`:

```python
import pandas as pd
from spectrum_annotator_ddzby import SpectrumAnnotator

# Load spectrum from CSV
spectrum_df = pd.read_csv("my_spectrum.csv")
exp_mz = spectrum_df['mz'].values
exp_intensity = spectrum_df['intensity'].values

# Create annotator and generate plot
annotator = SpectrumAnnotator(
    peptide="PEPTIDEK",
    modifications=[{'position': 5, 'residue': 'T', 'mass': 203.0794}],
    precursor_charge=2,
    precursor_mz=500.0,
    exp_mz=exp_mz,
    exp_intensity=exp_intensity,
    tolerance_ppm=20.0,
    site_index="TEST",
    gene="TEST"
)

fig = annotator.plot(output_path="annotated.pdf")
```

### Example 4: Batch Process Multiple Spectra

```python
from spectrum_annotator_ddzby import annotate_spectra_batch

# Process all spectra in a directory
output_files, stats_df = annotate_spectra_batch(
    summary_file="my_spectra_summary.csv",
    spectra_dir="spectra/",
    output_dir="output_pdfs/",
    tolerance_ppm=20.0,
    save_statistics=True
)

print(f"Generated {len(output_files)} annotated spectra")
print(f"Mean coverage: {stats_df['sequence_coverage'].mean()*100:.1f}%")
```

### Example 5: Access Glycan Library

```python
from spectrum_annotator_ddzby import (
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    MONOSACCHARIDE_MASSES,
)

# Print all O-glycans
print("O-Glycan Library:")
for name, glycan in O_GLYCAN_COMPOSITIONS.items():
    print(f"  {name}: {glycan.mass:.4f} Da")

# Get specific glycan mass
sialyl_t = O_GLYCAN_COMPOSITIONS['Sialyl-T']
print(f"\nSialyl-T mass: {sialyl_t.mass:.4f} Da")
print(f"Composition: {sialyl_t.composition}")

# Calculate custom glycan mass
custom_mass = (
    2 * MONOSACCHARIDE_MASSES['HexNAc'] +
    5 * MONOSACCHARIDE_MASSES['Hex'] +
    1 * MONOSACCHARIDE_MASSES['Fuc']
)
print(f"\nHexNAc2Hex5Fuc1 mass: {custom_mass:.4f} Da")
```

---

## Input and Output Formats

### Spectrum Data Format

**Option 1: Space/Tab-separated text**
```
204.087 1000
366.140 800
500.25 500
```

**Option 2: CSV file**
```csv
mz,intensity
204.087,1000
366.140,800
500.25,500
```

### Modification Format

Modifications are specified as a list of dictionaries:
```python
modifications = [
    {'position': 0, 'residue': 'N-term', 'mass': 229.1629},   # N-terminal
    {'position': 5, 'residue': 'S', 'mass': 203.0794},        # Residue 5
    {'position': -1, 'residue': 'C-term', 'mass': 100.0},     # C-terminal
]
```

- `position`: 0 = N-terminus, -1 = C-terminus, 1-indexed for residues
- `residue`: Amino acid letter or 'N-term'/'C-term'
- `mass`: Modification mass in Daltons

### Common Modification Masses

| Modification | Mass (Da) |
|-------------|-----------|
| TMT6plex | 229.1629 |
| O-GlcNAc (HexNAc) | 203.0794 |
| HexNAc + TMT | 528.2859 |
| Carbamidomethyl (Cys) | 57.0215 |
| Oxidation (Met) | 15.9949 |

### Output: Annotated Spectrum PDF

The PDF contains:
- Peptide sequence with fragmentation sites marked
- Sample information (gene, site, precursor m/z, charge)
- Sequence coverage and False Match Rate
- Color-coded annotated spectrum
- Mass error plot
- Legend

### Ion Color Scheme

| Ion | Color | Description |
|-----|-------|-------------|
| b | Blue | N-terminal HCD |
| y | Red | C-terminal HCD |
| c | Green | N-terminal ETD |
| z | Orange | C-terminal ETD |
| Y | Purple | Glycopeptide Y ions |
| Oxonium | Brown | Glycan diagnostic |
| Unassigned | Gray | No match |

---

## Troubleshooting

### "command not found: python" or "python is not recognized"

**Windows**: Python was not added to PATH during installation.
- Reinstall Python and check "Add Python to PATH"

**macOS/Linux**: Use `python3` instead of `python`

### "pip: command not found"

Use `pip3` instead of `pip`:
```bash
pip3 install spectrum-annotator-ddzby
```

### "Permission denied" when installing

Add `--user` flag:
```bash
pip3 install --user "spectrum-annotator-ddzby[streamlit]"
```

### Streamlit app won't start

1. Make sure you're in the correct directory:
   ```bash
   cd /path/to/GlycoSpectrumAnnotator
   ls app.py  # Should show app.py
   ```

2. Check if port 8501 is already in use:
   ```bash
   # macOS/Linux
   lsof -i :8501

   # Windows
   netstat -an | findstr 8501
   ```

3. Use a different port:
   ```bash
   streamlit run app.py --server.port 8502
   ```

### "ModuleNotFoundError: No module named 'spectrum_annotator_ddzby'"

The package wasn't installed correctly. Try:
```bash
pip3 install spectrum-annotator-ddzby
```

### No peaks are being matched

- Increase mass tolerance (try 30-50 ppm)
- Verify your modifications are correct
- Check that precursor charge is correct
- Make sure m/z values are in the expected range

### Memory error during batch processing

- Process fewer spectra at a time
- Close plots after saving: `plt.close(fig)`

---

## Reference

### Supported Glycans

**O-Glycans**: Tn, T-antigen, Core 1-4, Sialyl-Tn, Sialyl-T, Disialyl-T, and more

**N-Glycans**: Man5-Man9, A2G0, A2G1, A2G2, A2G2F, A2G2S1, A2G2S2, and more

### Supported Crosslinkers

| Name | Type | Stub Masses |
|------|------|-------------|
| DSSO | MS-cleavable | A: 54.01, T: 85.98, S: 103.99 |
| DSBSO | MS-cleavable | A: 54.01, T: 85.98, S: 103.99 |
| BS3 | Non-cleavable | - |
| DSS | Non-cleavable | - |

### Citation

If you use this tool in your research, please cite:

> Fu, L. (2026). Spectrum Annotator Ddzby: Universal MS/MS Spectrum Annotation Tool.
> GitHub: https://github.com/lfu46/GlycoSpectrumAnnotator

### Contact

- **Author**: Longping Fu
- **Email**: lpfu46@gatech.edu
- **GitHub Issues**: https://github.com/lfu46/GlycoSpectrumAnnotator/issues

### References

- Schulte, D. et al. "A Universal Spectrum Annotator for Complex Peptidoforms in Mass Spectrometry-Based Proteomics." *Analytical Chemistry* 2025, 97, 23120-23130.
