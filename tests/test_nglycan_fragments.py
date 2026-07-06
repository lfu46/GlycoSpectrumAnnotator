"""Regression tests for the N-glycan HCD fragment generators.

`generate_n_glycan_remainders` and `generate_n_glycan_bions` were promoted from
`nglyco_annotate_spectra` into the shared glycan_library (2026-07-06) so the
engine owns N-glycan fragment math alongside `generate_n_glycan_y_ions`. These
tests lock the composition-gating behavior that makes the ions trustworthy
(no NeuAc oxonium for asialo glycans, biosynthetically valid remainders).
"""
import pytest

from spectrum_annotator_ddzby import (
    generate_n_glycan_remainders,
    generate_n_glycan_bions,
)
from spectrum_annotator_ddzby.glycan_library import _classify_glycan_type_simple


HEXNAC = 203.0794
NEUAC = 291.0954


def test_classify_simple():
    assert _classify_glycan_type_simple({"HexNAc": 2, "Hex": 8}) == "high_mannose"
    assert _classify_glycan_type_simple({"HexNAc": 2, "Hex": 3}) == "paucimannose"
    assert _classify_glycan_type_simple({"HexNAc": 3, "Hex": 6}) == "hybrid"
    assert _classify_glycan_type_simple({"HexNAc": 4, "Hex": 5}) == "complex"
    assert _classify_glycan_type_simple({"HexNAc": 1, "Hex": 1}) == "unknown"


def test_remainders_have_y0_and_intact():
    rem = generate_n_glycan_remainders({"HexNAc": 2, "Hex": 8})
    assert rem[""] == 0.0                      # Y0 / bare peptide
    assert "intact" in rem                     # full glycan
    # reducing-end GlcNAc (Y1) at one HexNAc residue mass
    assert rem["+N"] == pytest.approx(HEXNAC, abs=1e-3)
    # no partial exceeds the intact mass
    assert max(rem.values()) == pytest.approx(rem["intact"], abs=1e-3)


def test_highmannose_mannose_ladder():
    rem = generate_n_glycan_remainders({"HexNAc": 2, "Hex": 8})
    # linear mannose tree: N2H4..N2H8 intermediates all present
    for h in range(4, 9):
        assert f"+N2H{h}" in rem


def test_bions_gate_on_sialic_acid():
    # asialo high-mannose must NOT produce a NeuAc oxonium
    asialo = generate_n_glycan_bions({"HexNAc": 2, "Hex": 8})
    assert "NeuAc" not in asialo
    # a sialylated complex glycan must produce NeuAc (292) and its -H2O (274)
    sialo = generate_n_glycan_bions({"HexNAc": 4, "Hex": 5, "NeuAc": 2})
    assert sialo["NeuAc"] == pytest.approx(NEUAC + 1.007276, abs=1e-3)
    assert "NeuAc-H2O" in sialo


def test_bions_lacdinac_gate():
    # LacdiNAc marker (HexNAc-HexNAc) only when extra HexNAc present (n_N >= 4)
    assert "HexNAc-HexNAc" in generate_n_glycan_bions({"HexNAc": 5, "Hex": 4})
    assert "HexNAc-HexNAc" not in generate_n_glycan_bions({"HexNAc": 3, "Hex": 5})


def test_bions_always_have_core():
    # every N-glycan carries the HexNAc + Hex core oxonium/disaccharide ions
    for comp in ({"HexNAc": 2, "Hex": 5}, {"HexNAc": 4, "Hex": 5, "NeuAc": 1}):
        ions = generate_n_glycan_bions(comp)
        assert "HexNAc" in ions and "Hex" in ions and "HexNAc-Hex" in ions
