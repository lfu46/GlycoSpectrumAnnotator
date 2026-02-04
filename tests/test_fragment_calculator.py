#!/usr/bin/env python3
"""
Tests for FragmentCalculator

Run with: pytest tests/test_fragment_calculator.py -v
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, '..')

from glycospectrum import (
    FragmentCalculator,
    TheoreticalIon,
    parse_modifications_from_string,
    match_peaks,
    AA_MASSES,
    MOD_MASSES,
)


class TestFragmentCalculator:
    """Tests for the FragmentCalculator class."""

    @pytest.fixture
    def simple_peptide(self):
        """Simple peptide without modifications."""
        return "PEPTIDE"

    @pytest.fixture
    def glycopeptide(self):
        """Glycopeptide with TMT and O-GlcNAc."""
        peptide = "AGYSQGATQYTQAQQTR"
        mods = parse_modifications_from_string("N-term(229.1629),4S(528.2859)")
        return peptide, mods

    def test_parse_modifications(self):
        """Test modification string parsing."""
        mod_string = "N-term(229.1629),4S(528.2859),17K(229.1629)"
        mods = parse_modifications_from_string(mod_string)

        assert len(mods) == 3
        assert mods[0]['position'] == 0  # N-term
        assert mods[1]['position'] == 4
        assert mods[1]['residue'] == 'S'
        assert abs(mods[1]['mass'] - 528.2859) < 0.001

    def test_precursor_mass(self, simple_peptide):
        """Test precursor mass calculation."""
        calc = FragmentCalculator(simple_peptide, [], 2)

        # Calculate expected mass manually
        expected = sum(AA_MASSES[aa] for aa in simple_peptide) + 18.010565  # + H2O

        assert abs(calc.precursor_mass - expected) < 0.01

    def test_b_ion_count(self, simple_peptide):
        """Test that correct number of b ions are generated."""
        calc = FragmentCalculator(simple_peptide, [], 2)
        b_ions = calc.calculate_b_ions(charges=[1])

        # Should have n-1 b ions for peptide of length n
        assert len(b_ions) == len(simple_peptide) - 1

    def test_y_ion_count(self, simple_peptide):
        """Test that correct number of y ions are generated."""
        calc = FragmentCalculator(simple_peptide, [], 2)
        y_ions = calc.calculate_y_ions(charges=[1])

        # Should have n-1 y ions for peptide of length n
        assert len(y_ions) == len(simple_peptide) - 1

    def test_glycopeptide_Y_ions(self, glycopeptide):
        """Test Y ion generation for glycopeptide."""
        peptide, mods = glycopeptide
        calc = FragmentCalculator(peptide, mods, 3)
        Y_ions = calc.calculate_Y_ions()

        # Should have Y0 (glycan loss) and Y1 (intact) at multiple charges
        assert len(Y_ions) > 0

        # Check for Y0 ion
        y0_ions = [ion for ion in Y_ions if ion.ion_number == 0]
        assert len(y0_ions) > 0

    def test_oxonium_ions(self, glycopeptide):
        """Test oxonium ion generation."""
        peptide, mods = glycopeptide
        calc = FragmentCalculator(peptide, mods, 3)
        oxonium = calc.calculate_oxonium_ions()

        # Should have multiple oxonium ions
        assert len(oxonium) > 0

        # All oxonium ions should be charge 1
        for ion in oxonium:
            assert ion.charge == 1

    def test_neutral_losses(self, simple_peptide):
        """Test neutral loss ion generation."""
        calc = FragmentCalculator(simple_peptide, [], 2)
        b_ions = calc.calculate_b_ions(charges=[1])
        nl_ions = calc.calculate_neutral_loss_ions(b_ions, ['H2O', 'NH3'])

        # Should have neutral loss variants
        assert len(nl_ions) > 0

        # Check that neutral loss is recorded
        for ion in nl_ions:
            assert ion.neutral_loss in ['H2O', 'NH3']


class TestPeakMatching:
    """Tests for peak matching functionality."""

    def test_exact_match(self):
        """Test matching with exact m/z values."""
        theo_ions = [
            TheoreticalIon('b', 1, 1, 100.0, 'P', '', 'b1'),
            TheoreticalIon('b', 2, 1, 200.0, 'PE', '', 'b2'),
        ]
        exp_mz = np.array([100.0, 200.0, 300.0])
        exp_intensity = np.array([1000.0, 2000.0, 500.0])

        matched = match_peaks(theo_ions, exp_mz, exp_intensity, tolerance_ppm=20.0)

        assert len(matched) == 2

    def test_tolerance(self):
        """Test that tolerance is correctly applied."""
        theo_ions = [
            TheoreticalIon('b', 1, 1, 100.0, 'P', '', 'b1'),
        ]

        # Within 20 ppm (100 * 20 / 1e6 = 0.002)
        exp_mz_close = np.array([100.001])
        exp_intensity = np.array([1000.0])

        matched = match_peaks(theo_ions, exp_mz_close, exp_intensity, tolerance_ppm=20.0)
        assert len(matched) == 1

        # Outside 20 ppm
        exp_mz_far = np.array([100.01])
        matched = match_peaks(theo_ions, exp_mz_far, exp_intensity, tolerance_ppm=20.0)
        assert len(matched) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
