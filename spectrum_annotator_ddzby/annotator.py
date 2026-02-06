#!/usr/bin/env python3
"""
Spectrum Annotator for Glycopeptides

This module creates publication-quality annotated EThcD spectra for
glycopeptide identification, similar to IPSA output.

Author: Longping Fu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from typing import List, Dict, Optional, Tuple, Union
import os
import json

# Import the fragment calculator (relative import for package)
from .fragment_calculator import (
    FragmentCalculator,
    TheoreticalIon,
    MatchedIon,
    FalseMatchRate,
    match_peaks,
    parse_modifications_from_string,
    calculate_false_match_rate,
    calculate_annotation_statistics,
    OXONIUM_IONS
)

# =============================================================================
# COLOR SCHEME (IPSA-style colors)
# =============================================================================

# Ion type colors (IPSA palette)
ION_COLORS = {
    'b': '#0d75bc',      # Blue (b ions - HCD N-terminal)
    'y': '#be202d',      # Red (y ions - HCD C-terminal)
    'c': '#07a14a',      # Green (c ions - ETD N-terminal)
    'z': '#f79420',      # Orange (z ions - ETD C-terminal)
    'Y': '#8491B4',      # Purple-gray (Y ions - intact glycopeptide)
    'Y0': '#9B59B6',     # Purple (Y0 - peptide with glycan loss)
    'Y1': '#8491B4',     # Purple-gray (Y1 - intact glycopeptide)
    'oxonium': '#7E6148', # Brown (oxonium/glycan diagnostic)
    'precursor': '#666666',  # Dark gray (precursor)
    'unassigned': '#a6a6a6', # Gray (unmatched - IPSA style)
}

# Neutral loss suffix colors (slightly lighter/transparent versions)
ION_COLORS_NL = {
    'b_NL': '#0d75bc80',
    'y_NL': '#be202d80',
    'c_NL': '#07a14a80',
    'z_NL': '#f7942080',
}

# =============================================================================
# SPECTRUM ANNOTATOR CLASS
# =============================================================================

class SpectrumAnnotator:
    """
    Creates annotated spectrum plots for glycopeptides.
    """

    def __init__(self,
                 peptide: str,
                 modifications: List[Dict],
                 precursor_charge: int,
                 precursor_mz: float,
                 exp_mz: np.ndarray,
                 exp_intensity: np.ndarray,
                 tolerance_ppm: float = 20.0,
                 site_index: str = "",
                 gene: str = "",
                 activation_type: str = "HCD"):
        """
        Initialize the spectrum annotator.

        Args:
            peptide: Peptide sequence
            modifications: List of modification dicts
            precursor_charge: Precursor charge state
            precursor_mz: Experimental precursor m/z
            exp_mz: Experimental m/z array
            exp_intensity: Experimental intensity array
            tolerance_ppm: Mass tolerance for matching (default 20 ppm)
            site_index: Site identifier (e.g., "Q96KR1_S195")
            gene: Gene name
            activation_type: Fragmentation method - "HCD" (b/y ions) or "EThcD" (b/y/c/z ions)
        """
        self.peptide = peptide
        self.modifications = modifications
        self.precursor_charge = precursor_charge
        self.precursor_mz = precursor_mz
        self.exp_mz = exp_mz
        self.exp_intensity = exp_intensity
        self.tolerance_ppm = tolerance_ppm
        self.site_index = site_index
        self.gene = gene
        self.activation_type = activation_type.upper()

        # Calculate theoretical fragments
        self.calculator = FragmentCalculator(
            peptide, modifications, precursor_charge, max_fragment_charge=2
        )

        # Get theoretical ions based on activation type
        if self.activation_type == "HCD":
            # HCD: only b, y, Y, precursor, oxonium (no c/z ions)
            # Note: NO glycan neutral losses - they don't provide localization info
            self.theoretical_ions = self.calculator.get_hcd_ions_flat(
                include_neutral_losses=True,
                neutral_loss_types=['H2O', 'NH3']  # No glycan losses
            )
        else:
            # EThcD or other: include all ions (b, y, c, z, Y, precursor, oxonium)
            # Note: NO glycan neutral losses - they don't provide localization info
            self.theoretical_ions = self.calculator.get_all_ions_flat(
                include_neutral_losses=True,
                neutral_loss_types=['H2O', 'NH3']  # No glycan losses
            )

        # Match peaks
        self.matched_ions = match_peaks(
            self.theoretical_ions, exp_mz, exp_intensity, tolerance_ppm
        )

        # Create lookup for matched peaks
        self.matched_mz_set = set()
        self.peak_annotations = {}  # exp_mz -> MatchedIon
        for ion in self.matched_ions:
            self.matched_mz_set.add(ion.exp_mz)
            # Keep the best annotation (highest intensity if multiple)
            if ion.exp_mz not in self.peak_annotations:
                self.peak_annotations[ion.exp_mz] = ion

        # Find glycan modification position
        # HexNAc+TMT = 528.2859 Da, HexNAc = 299.123 Da (metabolic labeling)
        self.glycan_position = None
        for mod in modifications:
            if abs(mod['mass'] - 528.2859) < 0.1 or abs(mod['mass'] - 299.123) < 0.1:
                self.glycan_position = mod['position']
                break

        # Calculate false match rate (spectrum shifting method from Schulte et al.)
        self.false_match_rate = calculate_false_match_rate(
            self.theoretical_ions, exp_mz, exp_intensity, tolerance_ppm
        )

        # Calculate comprehensive annotation statistics
        self.annotation_stats = calculate_annotation_statistics(
            self.matched_ions, self.theoretical_ions,
            exp_mz, exp_intensity, len(peptide)
        )

    def _get_ion_color(self, ion: 'MatchedIon', has_neutral_loss: bool = False) -> str:
        """Get color for an ion type."""
        ion_type = ion.ion_type

        # Special handling for Y ions based on ion_number (Y0 vs Y1)
        if ion_type == 'Y':
            if ion.ion_number == 0:  # Y0 - glycan loss
                return ION_COLORS.get('Y0', ION_COLORS['Y'])
            else:  # Y1 - intact glycopeptide
                return ION_COLORS.get('Y1', ION_COLORS['Y'])

        # Site-determining ions (with glycan) always use solid colors
        if ion.has_modification:
            return ION_COLORS.get(ion_type, ION_COLORS['unassigned'])

        if has_neutral_loss:
            return ION_COLORS_NL.get(f'{ion_type}_NL', ION_COLORS.get(ion_type, ION_COLORS['unassigned']))
        return ION_COLORS.get(ion_type, ION_COLORS['unassigned'])

    def _format_annotation(self, ion: MatchedIon, short: bool = False) -> str:
        """Format ion annotation for display."""
        if ion.ion_type == 'oxonium':
            return ion.annotation
        elif ion.ion_type == 'precursor':
            # Format as [M+nH]+n for precursor ions
            if 'charge_reduced' in ion.annotation.lower() or 'cr' in ion.annotation.lower():
                return f"[M+{ion.charge}H]+{ion.charge}" if short else ion.annotation
            return f"[M+{ion.charge}H]+{ion.charge}" if short else ion.annotation
        elif ion.ion_type == 'Y':
            # Y0 = peptide only (glycan lost), Y1 = intact glycopeptide
            y_type = "Y0" if ion.ion_number == 0 else "Y1"
            return f"{y_type} {ion.charge}+" if short else ion.annotation
        else:
            # b, y, c, z ions
            charge_str = "" if ion.charge == 1 else f" {ion.charge}+"  # Add space before charge
            # Add glycan marker if ion contains glycan (use short form for display)
            glycan_str = "+HexNAc" if ion.has_modification else ""
            nl_str = ""
            if ion.neutral_loss:
                nl_map = {'H2O': '°', 'NH3': '*', 'CO2': '**', 'HexNAc_TMT': '-Glyc', 'HexNAc': '-Hex'}
                nl_str = nl_map.get(ion.neutral_loss, f'-{ion.neutral_loss}')
            return f"{ion.ion_type}{ion.ion_number}{glycan_str}{nl_str}{charge_str}"

    def _get_fragmentation_coverage(self) -> Dict[str, set]:
        """Get which peptide bonds were fragmented by each ion type.

        Includes both base ions and neutral loss ions, since neutral loss
        ions also indicate bond cleavage.
        """
        coverage = {
            'b': set(),  # b ions (N-terminal, HCD)
            'c': set(),  # c ions (N-terminal, ETD)
            'y': set(),  # y ions (C-terminal, HCD)
            'z': set(),  # z ions (C-terminal, ETD)
        }

        for ion in self.matched_ions:
            if ion.ion_type in coverage:
                coverage[ion.ion_type].add(ion.ion_number)

        return coverage

    def _get_glycan_containing_ions(self) -> Dict[str, set]:
        """Get which matched ions contain the glycan (site-determining ions).

        Returns a dict mapping ion type to set of ion numbers that contain glycan.
        """
        glycan_ions = {
            'b': set(),
            'c': set(),
            'y': set(),
            'z': set(),
        }

        for ion in self.matched_ions:
            if ion.ion_type in glycan_ions and ion.has_modification:
                glycan_ions[ion.ion_type].add(ion.ion_number)

        return glycan_ions

    def plot(self,
             figsize: Tuple[float, float] = (7, 5),
             output_path: Optional[str] = None,
             show_error_plot: bool = True,
             intensity_threshold_pct: float = 1.0,
             max_labels: int = 50) -> plt.Figure:
        """
        Create annotated spectrum plot.

        Args:
            figsize: Figure size (width, height)
            output_path: Path to save PDF (optional)
            show_error_plot: Whether to show mass error plot
            intensity_threshold_pct: Minimum intensity (% base peak) to label
            max_labels: Maximum number of labels to show

        Returns:
            matplotlib Figure object
        """
        # Set Arial as default font for all text
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.sans-serif'] = ['Arial']

        # Normalize intensities to % base peak
        base_peak = np.max(self.exp_intensity)
        rel_intensity = self.exp_intensity / base_peak * 100

        # Set up figure
        if show_error_plot:
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(4, 1, height_ratios=[1.2, 0.5, 3, 1], hspace=0.05)
            ax_seq = fig.add_subplot(gs[0])
            ax_info = fig.add_subplot(gs[1])
            ax_spec = fig.add_subplot(gs[2])
            ax_error = fig.add_subplot(gs[3], sharex=ax_spec)
        else:
            fig = plt.figure(figsize=(figsize[0], figsize[1] * 0.8))
            gs = GridSpec(3, 1, height_ratios=[1.2, 0.5, 3], hspace=0.05)
            ax_seq = fig.add_subplot(gs[0])
            ax_info = fig.add_subplot(gs[1])
            ax_spec = fig.add_subplot(gs[2])
            ax_error = None

        # =====================================================================
        # 1. Peptide Sequence with Fragmentation Sites (IPSA-style)
        # =====================================================================
        ax_seq.set_xlim(-0.5, len(self.peptide) + 0.5)
        ax_seq.set_ylim(-0.5, 1.5)  # Expanded for ion labels
        ax_seq.axis('off')

        # Get fragmentation coverage by ion type
        coverage = self._get_fragmentation_coverage()
        # Get glycan-containing (site-determining) ions
        glycan_ions = self._get_glycan_containing_ions()

        # Draw peptide sequence (no boxes, Arial font, size 8)
        for i, aa in enumerate(self.peptide):
            # Highlight glycosylated residue with color only (no box)
            if self.glycan_position and i + 1 == self.glycan_position:
                text_color = '#D2691E'  # Chocolate color for glycosylated residue
                fontweight = 'bold'
            else:
                text_color = 'black'
                fontweight = 'normal'

            ax_seq.text(i + 0.5, 0.5, aa, ha='center', va='center', fontsize=8,
                       fontweight=fontweight, fontfamily='Arial', color=text_color)

            # Draw IPSA-style fragmentation marks (vertical + diagonal)
            # b/y have longer vertical, c/z have shorter vertical, same x position
            if i < len(self.peptide) - 1:
                bond_pos = i + 1  # N-terminal ion number (b1, c1 after first AA)
                c_bond = len(self.peptide) - bond_pos  # C-terminal ion number

                # Base position for marks (between amino acids)
                x_base = i + 0.85

                # N-terminal ions (b/c) - vertical from middle up, then diagonal up-LEFT
                # (b ions cover residues TO THE LEFT of cleavage)
                # b ions (blue) - longer vertical
                if bond_pos in coverage['b']:
                    y_top_b = 0.9  # longer vertical
                    # Use thicker line if glycan-containing (site-determining)
                    lw = 2.5 if bond_pos in glycan_ions['b'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_top_b], color=ION_COLORS['b'], linewidth=lw)  # vertical
                    ax_seq.plot([x_base, x_base - 0.15], [y_top_b, y_top_b + 0.10], color=ION_COLORS['b'], linewidth=lw)  # diagonal LEFT
                    # Add label for site-determining ions (b ions cover residues to the LEFT)
                    if bond_pos in glycan_ions['b']:
                        ax_seq.text(x_base - 0.6, y_top_b + 0.25, f'←b{bond_pos}+G', fontsize=5,
                                   color=ION_COLORS['b'], fontweight='bold', fontfamily='Arial')
                # c ions (green) - shorter vertical
                if bond_pos in coverage['c']:
                    y_top_c = 0.75  # shorter vertical
                    lw = 2.5 if bond_pos in glycan_ions['c'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_top_c], color=ION_COLORS['c'], linewidth=lw)  # vertical
                    ax_seq.plot([x_base, x_base - 0.15], [y_top_c, y_top_c + 0.10], color=ION_COLORS['c'], linewidth=lw)  # diagonal LEFT
                    if bond_pos in glycan_ions['c']:
                        ax_seq.text(x_base - 0.6, y_top_c + 0.15, f'←c{bond_pos}+G', fontsize=5,
                                   color=ION_COLORS['c'], fontweight='bold', fontfamily='Arial')

                # C-terminal ions (y/z) - vertical from middle down, then diagonal down-RIGHT
                # (y ions cover residues TO THE RIGHT of cleavage)
                # y ions (red) - longer vertical
                if c_bond in coverage['y']:
                    y_bot_y = 0.1  # longer vertical (goes lower)
                    lw = 2.5 if c_bond in glycan_ions['y'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_bot_y], color=ION_COLORS['y'], linewidth=lw)  # vertical
                    ax_seq.plot([x_base, x_base + 0.15], [y_bot_y, y_bot_y - 0.10], color=ION_COLORS['y'], linewidth=lw)  # diagonal RIGHT
                    if c_bond in glycan_ions['y']:
                        # For y ions, show label with arrow pointing right (→) to indicate
                        # the ion covers residues TO THE RIGHT of the cleavage site
                        ax_seq.text(x_base + 0.2, y_bot_y - 0.30, f'y{c_bond}+G→', fontsize=5,
                                   color=ION_COLORS['y'], fontweight='bold', fontfamily='Arial')
                # z ions (orange) - shorter vertical
                if c_bond in coverage['z']:
                    y_bot_z = 0.25  # shorter vertical (doesn't go as low)
                    lw = 2.5 if c_bond in glycan_ions['z'] else 1.5
                    ax_seq.plot([x_base, x_base], [0.5, y_bot_z], color=ION_COLORS['z'], linewidth=lw)  # vertical
                    ax_seq.plot([x_base, x_base + 0.15], [y_bot_z, y_bot_z - 0.10], color=ION_COLORS['z'], linewidth=lw)  # diagonal RIGHT
                    if c_bond in glycan_ions['z']:
                        ax_seq.text(x_base + 0.2, y_bot_z - 0.20, f'z{c_bond}+G→', fontsize=5,
                                   color=ION_COLORS['z'], fontweight='bold', fontfamily='Arial')

        # =====================================================================
        # 2. Info Panel (with False Match Rate)
        # =====================================================================
        ax_info.axis('off')

        # Calculate statistics
        total_bonds = len(self.peptide) - 1
        # Combine all N-terminal and C-terminal ion coverages
        n_bonds = coverage['b'] | coverage['c']
        c_bonds = coverage['y'] | coverage['z']
        fragmented_bonds = len(n_bonds | c_bonds)
        coverage_pct = fragmented_bonds / total_bonds * 100 if total_bonds > 0 else 0

        # Get false match rate values
        fmr = self.false_match_rate
        stats = self.annotation_stats

        # Line 1: Basic info
        info_line1 = (f"Gene: {self.gene}    Site: {self.site_index}    "
                     f"Precursor m/z: {self.precursor_mz:.4f}    "
                     f"Charge: +{self.precursor_charge}")

        # Line 2: Annotation statistics with FMR
        # Seq. Coverage = fraction of peptide bonds with detected fragment ions
        # FMR = False Match Rate (estimated spurious matches via spectrum shifting)
        info_line2 = (f"Seq. Coverage: {fragmented_bonds}/{total_bonds} ({coverage_pct:.0f}%)    "
                     f"FMR: {fmr.fmr_peaks*100:.1f}%")

        ax_info.text(0.5, 0.75, info_line1, ha='center', va='center', fontsize=9,
                    fontfamily='Arial', transform=ax_info.transAxes)
        ax_info.text(0.5, 0.25, info_line2, ha='center', va='center', fontsize=8,
                    fontfamily='Arial', transform=ax_info.transAxes, color='#555555')

        # =====================================================================
        # 3. Spectrum Plot
        # =====================================================================
        # Plot all peaks first (unassigned color)
        for i, (mz, inten) in enumerate(zip(self.exp_mz, rel_intensity)):
            if mz not in self.matched_mz_set:
                ax_spec.vlines(mz, 0, inten, color=ION_COLORS['unassigned'], linewidth=0.8)

        # Plot matched peaks with colors
        labels_added = 0
        labeled_positions = []  # Track label positions to avoid overlap

        # Sort matched ions for labeling priority:
        # 1. Site-determining ions (has_modification=True) get highest priority
        # 2. Then by intensity (highest first)
        def label_priority(ion):
            # Priority: glycan-containing ions first (0), then others (1)
            # Within each group, sort by intensity (descending, so negate)
            priority = 0 if ion.has_modification else 1
            return (priority, -ion.exp_intensity)

        sorted_matched = sorted(self.matched_ions, key=label_priority)

        for ion in sorted_matched:
            color = self._get_ion_color(ion, bool(ion.neutral_loss))
            rel_int = ion.exp_intensity / base_peak * 100

            # Plot the peak
            ax_spec.vlines(ion.exp_mz, 0, rel_int, color=color, linewidth=1.5)

            # Add label if above threshold and not too many labels
            if rel_int >= intensity_threshold_pct and labels_added < max_labels:
                # Check for label overlap (reduced threshold from 30 to 15 m/z)
                can_label = True
                for pos in labeled_positions:
                    if abs(ion.exp_mz - pos[0]) < 15 and abs(rel_int - pos[1]) < 8:
                        can_label = False
                        break

                if can_label:
                    label = self._format_annotation(ion, short=True)
                    ax_spec.annotate(label, (ion.exp_mz, rel_int),
                                    textcoords="offset points", xytext=(0, 5),
                                    ha='center', va='bottom', fontsize=7,
                                    fontfamily='Arial', color=color, fontweight='bold',
                                    rotation=90 if len(label) > 4 else 0)
                    labeled_positions.append((ion.exp_mz, rel_int))
                    labels_added += 1

        # Set spectrum axes
        ax_spec.set_xlim(min(self.exp_mz) - 50, max(self.exp_mz) + 50)
        ax_spec.set_ylim(0, 110)
        ax_spec.set_ylabel('Relative Abundance (%)', fontsize=9, fontfamily='Arial')
        ax_spec.spines['top'].set_visible(False)
        ax_spec.spines['right'].set_visible(False)
        ax_spec.tick_params(labelsize=8)

        if not show_error_plot:
            ax_spec.set_xlabel('m/z', fontsize=9, fontfamily='Arial')

        # =====================================================================
        # 4. Mass Error Plot
        # =====================================================================
        if ax_error is not None:
            for ion in self.matched_ions:
                color = self._get_ion_color(ion, bool(ion.neutral_loss))
                ax_error.scatter(ion.exp_mz, ion.mass_error_ppm, color=color, s=20, alpha=0.7)

            ax_error.axhline(0, color='black', linewidth=0.5, linestyle='--')
            ax_error.axhline(self.tolerance_ppm, color='red', linewidth=0.5, linestyle=':')
            ax_error.axhline(-self.tolerance_ppm, color='red', linewidth=0.5, linestyle=':')

            ax_error.set_ylim(-self.tolerance_ppm * 1.5, self.tolerance_ppm * 1.5)
            ax_error.set_xlabel('m/z', fontsize=9, fontfamily='Arial')
            ax_error.set_ylabel('Error (ppm)', fontsize=9, fontfamily='Arial')
            ax_error.spines['top'].set_visible(False)
            ax_error.spines['right'].set_visible(False)
            ax_error.tick_params(labelsize=8)

        # =====================================================================
        # 5. Legend
        # =====================================================================
        legend_elements = [
            Line2D([0], [0], color=ION_COLORS['b'], linewidth=2, label='b'),
            Line2D([0], [0], color=ION_COLORS['y'], linewidth=2, label='y'),
            Line2D([0], [0], color=ION_COLORS['c'], linewidth=2, label='c'),
            Line2D([0], [0], color=ION_COLORS['z'], linewidth=2, label='z'),
            Line2D([0], [0], color=ION_COLORS['Y0'], linewidth=2, label='Y0'),
            Line2D([0], [0], color=ION_COLORS['Y1'], linewidth=2, label='Y1'),
            Line2D([0], [0], color=ION_COLORS['oxonium'], linewidth=2, label='Oxonium'),
        ]
        ax_spec.legend(handles=legend_elements, loc='upper right', fontsize=7, ncol=4,
                      prop={'family': 'Arial', 'size': 7}, handlelength=1.2, handletextpad=0.4,
                      columnspacing=0.8, borderpad=0.4)

        # Add annotation symbol legend (neutral losses and glycan)
        symbol_legend = "+HexNAc: with glycan    °: -H2O    *: -NH3    **: -CO2"
        ax_spec.text(0.02, 0.97, symbol_legend, transform=ax_spec.transAxes,
                    fontsize=6, fontfamily='Arial', color='#555555',
                    verticalalignment='top', horizontalalignment='left')

        # Title
        fig.suptitle(f'{self.gene} - {self.site_index}', fontsize=10, fontweight='bold',
                    fontfamily='Arial', y=0.98)

        plt.tight_layout()

        # Save if output path provided
        if output_path:
            plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')
            print(f"Saved: {output_path}")

        return fig


# =============================================================================
# BATCH ANNOTATION FUNCTION
# =============================================================================

def annotate_spectra_batch(summary_file: str,
                          spectra_dir: str,
                          output_dir: str,
                          n_spectra: int = None,
                          tolerance_ppm: float = 20.0,
                          save_statistics: bool = True) -> Tuple[List[str], Optional[pd.DataFrame]]:
    """
    Batch annotate multiple spectra.

    Args:
        summary_file: Path to extracted spectra summary CSV
        spectra_dir: Directory containing spectrum CSV files
        output_dir: Directory for output PDFs
        n_spectra: Number of spectra to process (None = all)
        tolerance_ppm: Mass tolerance in ppm
        save_statistics: Whether to save annotation statistics CSV

    Returns:
        Tuple of (list of output file paths, DataFrame with statistics)
    """
    # Load summary
    df = pd.read_csv(summary_file)

    if n_spectra:
        df = df.head(n_spectra)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    output_files = []
    statistics_list = []

    for idx, row in df.iterrows():
        print(f"\nProcessing {idx + 1}/{len(df)}: {row['site_index']}")

        try:
            # Load spectrum data
            spec_file = os.path.join(spectra_dir, row['spectrum_file'])
            spec_df = pd.read_csv(spec_file)
            exp_mz = spec_df['mz'].values
            exp_intensity = spec_df['intensity'].values

            # Parse modifications
            if 'modifications_json' in row and pd.notna(row['modifications_json']):
                modifications = json.loads(row['modifications_json'])
            else:
                modifications = parse_modifications_from_string(row.get('Assigned_Modifications', ''))

            # Create annotator
            annotator = SpectrumAnnotator(
                peptide=row['Peptide'],
                modifications=modifications,
                precursor_charge=int(row['Charge']),
                precursor_mz=float(row.get('Calibrated_Observed_MZ', row.get('mzML_precursor_mz', 0))),
                exp_mz=exp_mz,
                exp_intensity=exp_intensity,
                tolerance_ppm=tolerance_ppm,
                site_index=row['site_index'],
                gene=row['Gene']
            )

            # Generate output filename
            output_file = os.path.join(output_dir, f"{row['site_index']}_{row['scan_number']}.pdf")

            # Create plot
            fig = annotator.plot(output_path=output_file)
            plt.close(fig)

            output_files.append(output_file)

            # Collect statistics
            fmr = annotator.false_match_rate
            stats = annotator.annotation_stats
            statistics_list.append({
                'site_index': row['site_index'],
                'gene': row['Gene'],
                'peptide': row['Peptide'],
                'charge': int(row['Charge']),
                'scan_number': row.get('scan_number', ''),
                'sequence_coverage': stats['sequence_coverage'],
                'sequence_coverage_bonds': stats['sequence_coverage_bonds'],
                'peaks_annotated': stats['peaks_annotated'],
                'peaks_annotated_count': stats['peaks_annotated_count'],
                'intensity_annotated': stats['intensity_annotated'],
                'fragments_found': stats['fragments_found'],
                'fragments_found_count': stats['fragments_found_count'],
                'fmr_peaks': fmr.fmr_peaks,
                'fmr_intensity': fmr.fmr_intensity,
                'matched_peaks': fmr.matched_peaks,
                'matched_intensity': fmr.matched_intensity,
                'avg_random_peaks': fmr.avg_random_peaks,
                'avg_random_intensity': fmr.avg_random_intensity,
            })

        except Exception as e:
            print(f"  Error: {e}")
            continue

    # Create statistics DataFrame
    stats_df = pd.DataFrame(statistics_list) if statistics_list else None

    # Save statistics if requested
    if save_statistics and stats_df is not None:
        stats_file = os.path.join(output_dir, "annotation_statistics.csv")
        stats_df.to_csv(stats_file, index=False)
        print(f"\nSaved annotation statistics to: {stats_file}")

    return output_files, stats_df


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Test with 10 example spectra from HEK293T
    SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
    OUTPUT_PATH = os.path.join(SOURCE_PATH, "point_to_point_response/")

    summary_file = os.path.join(OUTPUT_PATH, "extracted_spectra_EThcD_HEK293T_summary.csv")
    spectra_dir = os.path.join(OUTPUT_PATH, "extracted_spectra_EThcD_HEK293T")
    output_dir = os.path.join(OUTPUT_PATH, "annotated_spectra_test")

    print("=" * 70)
    print("Spectrum Annotation Test - 10 Example Spectra")
    print("(with False Match Rate calculation)")
    print("=" * 70)

    output_files, stats_df = annotate_spectra_batch(
        summary_file=summary_file,
        spectra_dir=spectra_dir,
        output_dir=output_dir,
        n_spectra=10,
        tolerance_ppm=20.0,
        save_statistics=True
    )

    print("\n" + "=" * 70)
    print(f"Generated {len(output_files)} annotated spectra")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

    # Display summary statistics
    if stats_df is not None and len(stats_df) > 0:
        print("\nAnnotation Statistics Summary:")
        print("-" * 50)
        print(f"Mean Sequence Coverage: {stats_df['sequence_coverage'].mean()*100:.1f}%")
        print(f"Mean Intensity Annotated: {stats_df['intensity_annotated'].mean()*100:.1f}%")
        print(f"Mean FMR (peaks): {stats_df['fmr_peaks'].mean()*100:.1f}%")
        print(f"Mean FMR (intensity): {stats_df['fmr_intensity'].mean()*100:.1f}%")
        print("-" * 50)
        print("\nPer-spectrum statistics:")
        print(stats_df[['site_index', 'sequence_coverage', 'intensity_annotated',
                       'fmr_peaks', 'fmr_intensity']].to_string(index=False))
