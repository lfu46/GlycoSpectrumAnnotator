#!/usr/bin/env python3
"""
Extract EThcD spectrum data from calibrated mzML files for O-GlcNAc PSM annotation.
This script reads the EThcD ranked CSV files and extracts complete spectrum information
including m/z, intensity, precursor info, and metadata needed for spectral annotation.
"""

import os
import pandas as pd
import numpy as np
from pyteomics import mzml
from collections import defaultdict
import json

# Configuration
SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_PATH = os.path.join(SOURCE_PATH, "point_to_point_response/")

# mzML directories for each cell type
MZML_DIRS = {
    "HEK293T": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/",
    "HepG2": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/",
    "Jurkat": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/"
}

# EThcD ranked file paths
ETHCD_FILES = {
    "HEK293T": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_HEK293T_EThcD_ranked.csv"),
    "HepG2": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_HepG2_EThcD_ranked.csv"),
    "Jurkat": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_Jurkat_EThcD_ranked.csv")
}


def parse_spectrum_id(spectrum_str):
    """
    Parse spectrum identifier to extract file name and scan number.
    Format: {file_name}.{scan}.{scan}.{charge}
    """
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        file_name = parts[0]
        scan_number = int(parts[1])
        charge = int(parts[3])
        return file_name, scan_number, charge
    return None, None, None


def find_calibrated_mzml(file_name, mzml_dir):
    """Find the calibrated mzML file path for a given file name."""
    patterns = [
        f"{file_name}_calibrated.mzML",
        f"{file_name}_mz_calibrated.mzML",
        f"{file_name}_ppm_calibrated.mzML",
    ]

    for pattern in patterns:
        path = os.path.join(mzml_dir, pattern)
        if os.path.exists(path):
            return path

    # Search for any matching calibrated file
    try:
        all_files = os.listdir(mzml_dir)
        for f in all_files:
            if f.startswith(file_name) and 'calibrated' in f.lower() and f.endswith('.mzML'):
                return os.path.join(mzml_dir, f)
    except:
        pass

    return None


def extract_spectrum_data(mzml_reader, scan_number):
    """
    Extract complete spectrum data for a specific scan using indexed access.

    Returns dict with:
        - mz_array: numpy array of m/z values
        - intensity_array: numpy array of intensity values
        - precursor_mz: precursor m/z
        - precursor_charge: precursor charge state
        - precursor_intensity: precursor intensity
        - ms_level: MS level (should be 2)
        - retention_time: retention time in minutes
        - tic: total ion current
        - base_peak_mz: m/z of base peak
        - base_peak_intensity: intensity of base peak
        - filter_string: instrument filter string
    """
    scan_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"

    try:
        spectrum = mzml_reader.get_by_id(scan_id)
    except KeyError:
        # Try to find by scan number in ID
        for spec_id in mzml_reader.index.keys():
            if f"scan={scan_number}" in spec_id:
                spectrum = mzml_reader.get_by_id(spec_id)
                break
        else:
            return None

    result = {
        'mz_array': spectrum.get('m/z array', np.array([])),
        'intensity_array': spectrum.get('intensity array', np.array([])),
        'ms_level': spectrum.get('ms level'),
        'tic': spectrum.get('total ion current'),
        'base_peak_mz': spectrum.get('base peak m/z'),
        'base_peak_intensity': spectrum.get('base peak intensity'),
    }

    # Get retention time and filter string from scanList
    if 'scanList' in spectrum and 'scan' in spectrum['scanList']:
        scan_info = spectrum['scanList']['scan'][0]
        rt = scan_info.get('scan start time')
        if rt is not None:
            # Convert to minutes if in seconds
            result['retention_time'] = float(rt)
        result['filter_string'] = scan_info.get('filter string', '')

    # Get precursor info for MS2 spectra
    if 'precursorList' in spectrum:
        precursor = spectrum['precursorList']['precursor'][0]
        if 'selectedIonList' in precursor:
            sel_ion = precursor['selectedIonList']['selectedIon'][0]
            result['precursor_mz'] = sel_ion.get('selected ion m/z')
            result['precursor_charge'] = sel_ion.get('charge state')
            result['precursor_intensity'] = sel_ion.get('peak intensity')

    return result


def parse_modifications(assigned_mods):
    """
    Parse the Assigned.Modifications string to extract modification positions and masses.

    Example: "N-term(229.1629),4S(528.2859),19K(229.1629)"

    Returns: list of dicts with 'position', 'residue', 'mass'
             position = 0 for N-term, -1 for C-term, 1-based for residues
    """
    if pd.isna(assigned_mods) or not assigned_mods:
        return []

    mods = []
    for mod in assigned_mods.split(','):
        mod = mod.strip()
        if not mod:
            continue

        # Extract mass from parentheses
        if '(' in mod and ')' in mod:
            mass_str = mod[mod.find('(')+1:mod.find(')')]
            try:
                mass = float(mass_str)
            except ValueError:
                continue

            position_part = mod[:mod.find('(')]

            if position_part == 'N-term':
                mods.append({'position': 0, 'residue': 'N-term', 'mass': mass})
            elif position_part == 'C-term':
                mods.append({'position': -1, 'residue': 'C-term', 'mass': mass})
            else:
                # Parse position and residue (e.g., "4S" or "19K")
                pos = ''
                res = ''
                for char in position_part:
                    if char.isdigit():
                        pos += char
                    else:
                        res += char
                if pos and res:
                    mods.append({'position': int(pos), 'residue': res, 'mass': mass})

    return mods


def extract_ethcd_spectra(cell_type):
    """
    Extract all EThcD spectra for a given cell type.

    Creates:
    - Individual spectrum CSV files with m/z and intensity
    - A summary CSV with all metadata
    - A JSON file with modification details
    """
    print(f"\n{'='*60}")
    print(f"Processing {cell_type} EThcD spectra")
    print(f"{'='*60}")

    # Load EThcD ranked file
    ethcd_file = ETHCD_FILES[cell_type]
    if not os.path.exists(ethcd_file):
        print(f"EThcD file not found: {ethcd_file}")
        return None

    df = pd.read_csv(ethcd_file)
    print(f"Loaded {len(df)} EThcD PSMs")

    mzml_dir = MZML_DIRS[cell_type]

    # Create output directory
    spectra_output_dir = os.path.join(OUTPUT_PATH, f"extracted_spectra_EThcD_{cell_type}")
    os.makedirs(spectra_output_dir, exist_ok=True)
    print(f"Output directory: {spectra_output_dir}")

    # Group PSMs by mzML file
    file_groups = defaultdict(list)
    for idx, row in df.iterrows():
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[file_name].append((idx, scan_number, row))

    print(f"PSMs from {len(file_groups)} unique mzML files")

    # Results storage
    results = []
    extracted_count = 0

    # Process each mzML file
    for file_name, scans in file_groups.items():
        mzml_path = find_calibrated_mzml(file_name, mzml_dir)

        if mzml_path is None:
            print(f"  WARNING: Calibrated mzML not found for: {file_name}")
            continue

        print(f"  Processing: {os.path.basename(mzml_path)} ({len(scans)} scans)...", end=" ", flush=True)

        try:
            with mzml.MzML(mzml_path, use_index=True) as reader:
                scans_extracted = 0

                for idx, scan_number, row in scans:
                    spec_data = extract_spectrum_data(reader, scan_number)

                    if spec_data is None:
                        continue

                    # Create unique filename
                    spec_filename = f"{row['site_index']}_{scan_number}.csv"
                    spec_filepath = os.path.join(spectra_output_dir, spec_filename)

                    # Save spectrum data (m/z and intensity)
                    spec_df = pd.DataFrame({
                        'mz': spec_data['mz_array'],
                        'intensity': spec_data['intensity_array']
                    })
                    spec_df.to_csv(spec_filepath, index=False)

                    # Parse modifications
                    mods = parse_modifications(row.get('Assigned.Modifications', ''))

                    # Collect result info
                    result = {
                        # Identifiers
                        'EThcD_Rank': row['EThcD_Rank'],
                        'Spectrum': row['Spectrum'],
                        'site_index': row['site_index'],
                        'scan_number': scan_number,
                        'spectrum_file': spec_filename,

                        # Peptide info
                        'Peptide': row['Peptide'],
                        'Modified_Peptide': row.get('Modified Peptide', ''),
                        'Assigned_Modifications': row.get('Assigned.Modifications', ''),
                        'modifications_json': json.dumps(mods),

                        # Protein info
                        'Gene': row['Gene'],
                        'Protein_ID': row['Protein.ID'],

                        # Site info
                        'peptide_site': row['peptide_site'],
                        'modified_residue': row['modified_residue'],
                        'site_number': row['site_number'],

                        # Glycan info
                        'Total_Glycan_Composition': row.get('Total.Glycan.Composition', ''),

                        # Precursor info from CSV
                        'Charge': row['Charge'],
                        'Observed_MZ': row.get('Observed.M.Z', None),
                        'Calibrated_Observed_MZ': row.get('Calibrated Observed M/Z', None),
                        'Calculated_Peptide_Mass': row.get('Calculated.Peptide.Mass', None),

                        # Precursor info from mzML
                        'mzML_precursor_mz': spec_data.get('precursor_mz'),
                        'mzML_precursor_charge': spec_data.get('precursor_charge'),
                        'mzML_precursor_intensity': spec_data.get('precursor_intensity'),

                        # Spectrum metadata
                        'retention_time': spec_data.get('retention_time'),
                        'tic': spec_data.get('tic'),
                        'base_peak_mz': spec_data.get('base_peak_mz'),
                        'base_peak_intensity': spec_data.get('base_peak_intensity'),
                        'n_peaks': len(spec_data['mz_array']),
                        'ms_level': spec_data.get('ms_level'),
                        'filter_string': spec_data.get('filter_string', ''),

                        # Scoring info
                        'Composite_Score': row['Composite_Score'],
                        'O_Pair_Score': row.get('O.Pair.Score', None),
                        'Hyperscore': row.get('Hyperscore', None),
                        'Confidence_Level': row.get('Confidence.Level', ''),

                        # mzML file info
                        'mzml_file': os.path.basename(mzml_path)
                    }

                    results.append(result)
                    scans_extracted += 1

                print(f"extracted {scans_extracted}/{len(scans)}")
                extracted_count += scans_extracted

        except Exception as e:
            print(f"Error: {e}")

    # Save results summary
    results_df = pd.DataFrame(results)
    summary_file = os.path.join(OUTPUT_PATH, f"extracted_spectra_EThcD_{cell_type}_summary.csv")
    results_df.to_csv(summary_file, index=False)

    print(f"\nExtracted {extracted_count}/{len(df)} EThcD spectra")
    print(f"Summary saved to: {summary_file}")
    print(f"Spectra saved to: {spectra_output_dir}/")

    return results_df


if __name__ == "__main__":
    import sys

    # Process all cell types by default, or a specific one if provided
    if len(sys.argv) > 1:
        cell_types = [sys.argv[1]]
    else:
        cell_types = ["HEK293T", "HepG2", "Jurkat"]

    print("Extracting EThcD spectra for O-GlcNAc annotation")
    print("=" * 60)

    all_results = {}
    for cell_type in cell_types:
        if cell_type not in MZML_DIRS:
            print(f"Unknown cell type: {cell_type}")
            continue

        results = extract_ethcd_spectra(cell_type)
        if results is not None:
            all_results[cell_type] = results

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for cell_type, results in all_results.items():
        print(f"{cell_type}: {len(results)} spectra extracted")

    print("\nDone!")
