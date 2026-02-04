#!/usr/bin/env python3
"""
Example: Annotate a single glycopeptide spectrum

This script demonstrates how to use GlycoSpectrumAnnotator to annotate
a single MS/MS spectrum of a glycopeptide.
"""

import sys
sys.path.insert(0, '..')

import pandas as pd
import numpy as np
from glycospectrum import (
    SpectrumAnnotator,
    FragmentCalculator,
    parse_modifications_from_string,
    match_peaks,
    calculate_false_match_rate,
)


def main():
    # Example peptide with O-GlcNAc modification
    peptide = "AGYSQGATQYTQAQQTR"
    mod_string = "N-term(229.1629),4S(528.2859)"  # TMT + HexNAc-TMT
    precursor_charge = 3

    # Parse modifications
    modifications = parse_modifications_from_string(mod_string)

    print("=" * 60)
    print("GlycoSpectrumAnnotator Example")
    print("=" * 60)
    print(f"Peptide: {peptide}")
    print(f"Modifications: {mod_string}")
    print(f"Charge: +{precursor_charge}")
    print()

    # Calculate theoretical fragments
    calc = FragmentCalculator(peptide, modifications, precursor_charge)

    print(f"Precursor mass: {calc.precursor_mass:.4f} Da")
    print(f"Precursor m/z: {calc.precursor_mz:.4f}")
    print()

    # Get all theoretical ions
    all_ions = calc.calculate_all_ions()

    print("Theoretical ions:")
    for ion_type, ions in all_ions.items():
        if ions:
            print(f"  {ion_type}: {len(ions)} ions")

    # If you have experimental data, you can annotate it:
    #
    # spectrum = pd.read_csv("your_spectrum.csv")
    # exp_mz = spectrum['mz'].values
    # exp_intensity = spectrum['intensity'].values
    #
    # annotator = SpectrumAnnotator(
    #     peptide=peptide,
    #     modifications=modifications,
    #     precursor_charge=precursor_charge,
    #     precursor_mz=calc.precursor_mz,
    #     exp_mz=exp_mz,
    #     exp_intensity=exp_intensity,
    #     tolerance_ppm=20.0,
    #     gene="GENE1",
    #     site_index="GENE1_S4"
    # )
    #
    # fig = annotator.plot(output_path="annotated_spectrum.pdf")
    # print(f"\nFalse Match Rate: {annotator.false_match_rate.fmr_peaks*100:.1f}%")

    print()
    print("=" * 60)
    print("Example complete!")
    print("To annotate real spectra, provide experimental m/z and intensity data.")
    print("=" * 60)


if __name__ == "__main__":
    main()
