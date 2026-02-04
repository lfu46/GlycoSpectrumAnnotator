#!/usr/bin/env python3
"""
Tests for glycan_library module

Run with: pytest tests/test_glycan_library.py -v
"""

import pytest
import sys
sys.path.insert(0, '..')

from glycospectrum import (
    MONOSACCHARIDE_MASSES,
    OXONIUM_IONS_EXTENDED,
    OXONIUM_IONS_O_GLCNAC,
    OXONIUM_IONS_N_GLYCAN,
    GlycanComposition,
    O_GLYCAN_COMPOSITIONS,
    N_GLYCAN_COMPOSITIONS,
    generate_y_ion_series,
    generate_n_glycan_y_ions,
    identify_glycan_from_mass,
    get_glycan_mass,
    # Crosslinker support
    Crosslinker,
    CROSSLINKERS,
    DSSO,
    DSBSO,
    BS3,
    generate_crosslink_fragments,
    identify_crosslink_stubs,
)


class TestMonosaccharideMasses:
    """Tests for monosaccharide mass definitions."""

    def test_hexnac_mass(self):
        """Test HexNAc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['HexNAc'] - 203.0794) < 0.001

    def test_hex_mass(self):
        """Test Hex mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['Hex'] - 162.0528) < 0.001

    def test_neuac_mass(self):
        """Test NeuAc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['NeuAc'] - 291.0954) < 0.001

    def test_fuc_mass(self):
        """Test Fuc mass is correct."""
        assert abs(MONOSACCHARIDE_MASSES['Fuc'] - 146.0579) < 0.001


class TestGlycanComposition:
    """Tests for GlycanComposition class."""

    def test_from_string_simple(self):
        """Test parsing simple composition string."""
        glycan = GlycanComposition.from_string("HexNAc1")
        assert glycan.composition == {'HexNAc': 1}
        assert abs(glycan.mass - 203.0794) < 0.01

    def test_from_string_complex(self):
        """Test parsing complex composition string."""
        glycan = GlycanComposition.from_string("HexNAc2Hex5")
        assert glycan.composition == {'HexNAc': 2, 'Hex': 5}
        expected_mass = 2 * 203.0794 + 5 * 162.0528
        assert abs(glycan.mass - expected_mass) < 0.01

    def test_from_string_with_sialic(self):
        """Test parsing composition with sialic acid."""
        glycan = GlycanComposition.from_string("HexNAc4Hex5NeuAc2")
        assert glycan.composition == {'HexNAc': 4, 'Hex': 5, 'NeuAc': 2}


class TestOxoniumIons:
    """Tests for oxonium ion definitions."""

    def test_basic_oxonium_present(self):
        """Test that basic HexNAc oxonium ions are defined."""
        assert 'HexNAc' in OXONIUM_IONS_O_GLCNAC
        assert 'HexNAc-H2O' in OXONIUM_IONS_O_GLCNAC

    def test_extended_includes_neuac(self):
        """Test that extended set includes NeuAc ions."""
        assert 'NeuAc' in OXONIUM_IONS_EXTENDED
        assert 'NeuAc-H2O' in OXONIUM_IONS_EXTENDED

    def test_n_glycan_includes_disaccharides(self):
        """Test that N-glycan set includes disaccharide ions."""
        assert 'HexNAc-Hex' in OXONIUM_IONS_N_GLYCAN


class TestYIonGeneration:
    """Tests for Y ion series generation."""

    def test_o_glcnac_y_ions(self):
        """Test Y ion generation for O-GlcNAc."""
        glycan = O_GLYCAN_COMPOSITIONS['O-GlcNAc']
        peptide_mass = 1000.0

        y_ions = generate_y_ion_series(glycan, peptide_mass)

        assert 'Y0' in y_ions
        assert 'Y1' in y_ions
        assert abs(y_ions['Y0'] - peptide_mass) < 0.01

    def test_n_glycan_y_ions(self):
        """Test Y ion generation for N-glycan."""
        composition = {'HexNAc': 2, 'Hex': 5}
        peptide_mass = 1500.0

        y_ions = generate_n_glycan_y_ions(composition, peptide_mass)

        assert 'Y0' in y_ions
        assert 'Y1' in y_ions
        assert 'Y2' in y_ions
        assert 'Y(core)' in y_ions

    def test_y_ions_include_fucose(self):
        """Test that fucosylated variants are generated."""
        composition = {'HexNAc': 2, 'Hex': 5, 'Fuc': 1}
        peptide_mass = 1500.0

        y_ions = generate_n_glycan_y_ions(composition, peptide_mass)

        # Should have fucosylated variants
        fuc_variants = [k for k in y_ions.keys() if 'F' in k]
        assert len(fuc_variants) > 0


class TestGlycanIdentification:
    """Tests for glycan identification from mass."""

    def test_identify_o_glcnac(self):
        """Test identifying O-GlcNAc from mass."""
        matches = identify_glycan_from_mass(203.08, tolerance_da=0.1)
        names = [m[0] for m in matches]
        assert 'O-GlcNAc' in names or 'O-GalNAc' in names

    def test_identify_n_glycan(self):
        """Test identifying N-glycan from mass."""
        # Man5 mass ~1216 Da
        matches = identify_glycan_from_mass(1216.42, tolerance_da=0.1)
        names = [m[0] for m in matches]
        assert 'Man5' in names


class TestGetGlycanMass:
    """Tests for get_glycan_mass function."""

    def test_simple_composition(self):
        """Test mass calculation from composition string."""
        mass = get_glycan_mass("HexNAc1")
        assert abs(mass - 203.0794) < 0.01

    def test_complex_composition(self):
        """Test mass calculation for complex composition."""
        mass = get_glycan_mass("HexNAc2Hex3")
        expected = 2 * 203.0794 + 3 * 162.0528
        assert abs(mass - expected) < 0.01


class TestExpandedOGlycans:
    """Tests for expanded O-glycan database (MSFragger 12 common structures)."""

    def test_all_12_msfragger_glycans_present(self):
        """Test that all 12 MSFragger O-glycans are defined."""
        # MSFragger glycan compositions
        expected_compositions = [
            {'HexNAc': 1},  # 1
            {'HexNAc': 1, 'Hex': 1},  # 2
            {'HexNAc': 1, 'NeuAc': 1},  # 3
            {'HexNAc': 2, 'Hex': 1},  # 4
            {'HexNAc': 1, 'Hex': 1, 'NeuAc': 1},  # 5
            {'HexNAc': 2, 'Hex': 2},  # 6
            {'HexNAc': 2, 'Hex': 1, 'NeuAc': 1},  # 7
            {'HexNAc': 1, 'Hex': 1, 'NeuAc': 2},  # 8
            {'HexNAc': 2, 'Hex': 2, 'NeuAc': 1},  # 9
            {'HexNAc': 2, 'Hex': 2, 'Fuc': 1, 'NeuAc': 1},  # 10
            {'HexNAc': 2, 'Hex': 2, 'NeuAc': 2},  # 11
            {'HexNAc': 2, 'Hex': 2, 'Fuc': 1, 'NeuAc': 2},  # 12
        ]

        for expected_comp in expected_compositions:
            found = False
            for glycan in O_GLYCAN_COMPOSITIONS.values():
                if glycan.composition == expected_comp:
                    found = True
                    break
            assert found, f"Missing O-glycan composition: {expected_comp}"

    def test_sialyl_tn_mass(self):
        """Test Sialyl-Tn (HexNAc1NeuAc1) mass."""
        glycan = O_GLYCAN_COMPOSITIONS['Sialyl-Tn']
        expected = MONOSACCHARIDE_MASSES['HexNAc'] + MONOSACCHARIDE_MASSES['NeuAc']
        assert abs(glycan.mass - expected) < 0.01

    def test_disialyl_t_mass(self):
        """Test Disialyl-T (HexNAc1Hex1NeuAc2) mass."""
        glycan = O_GLYCAN_COMPOSITIONS['Disialyl-T']
        expected = (MONOSACCHARIDE_MASSES['HexNAc'] +
                   MONOSACCHARIDE_MASSES['Hex'] +
                   2 * MONOSACCHARIDE_MASSES['NeuAc'])
        assert abs(glycan.mass - expected) < 0.01

    def test_core_structures_present(self):
        """Test that Core 1-4 structures are defined."""
        assert 'Core1' in O_GLYCAN_COMPOSITIONS
        assert 'Core2' in O_GLYCAN_COMPOSITIONS
        assert 'Core3' in O_GLYCAN_COMPOSITIONS
        assert 'Core4' in O_GLYCAN_COMPOSITIONS

    def test_phosphorylated_glycans(self):
        """Test phosphorylated O-glycan definitions."""
        assert 'HexNAc1-Phosphate' in O_GLYCAN_COMPOSITIONS
        glycan = O_GLYCAN_COMPOSITIONS['HexNAc1-Phosphate']
        assert 'Phosphate' in glycan.composition

    def test_sulfated_glycans(self):
        """Test sulfated O-glycan definitions."""
        assert 'HexNAc1-Sulfate' in O_GLYCAN_COMPOSITIONS
        glycan = O_GLYCAN_COMPOSITIONS['HexNAc1-Sulfate']
        assert 'Sulfate' in glycan.composition


class TestCrosslinkers:
    """Tests for crosslinker support."""

    def test_crosslinker_definitions(self):
        """Test that main crosslinkers are defined."""
        assert 'DSSO' in CROSSLINKERS
        assert 'DSBSO' in CROSSLINKERS
        assert 'BS3' in CROSSLINKERS

    def test_dsso_properties(self):
        """Test DSSO crosslinker properties."""
        assert DSSO.cleavable == True
        assert 'alkene' in DSSO.stub_masses
        assert 'thiol' in DSSO.stub_masses
        assert abs(DSSO.stub_masses['alkene'] - 54.0106) < 0.001

    def test_dsbso_properties(self):
        """Test DSBSO crosslinker properties."""
        assert DSBSO.cleavable == True
        assert DSBSO.reactive_groups == 'NHS'
        assert len(DSBSO.stub_masses) > 0

    def test_bs3_non_cleavable(self):
        """Test BS3 is non-cleavable."""
        assert BS3.cleavable == False
        assert len(BS3.stub_masses) == 0

    def test_stub_mass_difference(self):
        """Test A-T mass difference is ~32 Da (sulfur)."""
        alkene = DSSO.stub_masses['alkene']
        thiol = DSSO.stub_masses['thiol']
        diff = thiol - alkene
        assert abs(diff - 31.97) < 0.01


class TestCrosslinkFragments:
    """Tests for crosslink fragment generation."""

    def test_generate_dsso_fragments(self):
        """Test DSSO fragment generation."""
        pep1_mass = 1000.0
        pep2_mass = 1200.0

        fragments = generate_crosslink_fragments(
            pep1_mass, pep2_mass, DSSO, precursor_charge=4
        )

        assert 'peptide1' in fragments
        assert 'peptide2' in fragments
        assert 'α-A' in fragments['peptide1']
        assert 'β-T' in fragments['peptide2']

    def test_non_cleavable_fragments(self):
        """Test BS3 (non-cleavable) fragment generation."""
        pep1_mass = 1000.0
        pep2_mass = 1200.0

        fragments = generate_crosslink_fragments(
            pep1_mass, pep2_mass, BS3, precursor_charge=3
        )

        assert 'precursor' in fragments
        assert 'intact' in fragments['precursor']

    def test_identify_stubs(self):
        """Test stub identification from masses."""
        peptide_mass = 1000.0
        # Simulate observed masses with alkene and thiol stubs
        observed = [
            1000.0 + 54.01,  # Alkene stub
            1000.0 + 85.98,  # Thiol stub
            500.0,  # Random peak
        ]

        matches = identify_crosslink_stubs(
            observed, peptide_mass, DSSO, tolerance_da=0.05
        )

        assert len(matches) >= 2
        stub_names = [m[0] for m in matches]
        assert 'alkene' in stub_names
        assert 'thiol' in stub_names


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
