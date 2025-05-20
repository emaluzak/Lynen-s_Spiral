import pytest
from rdkit import Chem
import sys
import os

sys.path.append(os.path.abspath(r"src\lynen_spiral\Lynen_spiral_visualisation"))
from enhanced_fatty_acid import EnhancedFattyAcidMetabolism

# Type "pytest tests\test_enhanced_fatty_acid.py -v" to run tests with detailed output

# A set of fatty acids and edge cases
TEST_MOLECULES = {
    "saturated_even": "CCCCCCCCCC(=O)O",
    "unsaturated_trans": "C/C=C/CCCCCCCCC(=O)O",
    "unsaturated_cis": "C/C=C\\CCCCCCCC(=O)O",
    "odd_chain": "CCCCCCC(=O)O",
    "polyunsaturated": "CCC=CCC=CCC=CCC(=O)O",
    "short_chain": "CCCC(=O)O"
}
@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_initialization(name, smiles):
    if name == "not_fatty_acid":
        with pytest.raises(ValueError):
            EnhancedFattyAcidMetabolism(smiles)
    else:
        efa = EnhancedFattyAcidMetabolism(smiles)
        assert isinstance(efa.molecule, Chem.Mol)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_activation_and_transport(name, smiles):
    if name == "not_fatty_acid":
        return
    efa = EnhancedFattyAcidMetabolism(smiles)
    activated = efa.activate_fatty_acid()
    assert isinstance(activated, Chem.Mol)
    transported = efa.carnitine_shuttle(activated)
    assert isinstance(transported, Chem.Mol)
    assert "Transport" in efa.reaction_steps

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_alpha_beta_gamma_identification(name, smiles):
    if name == "not_fatty_acid" or name == "no_carbonyl":
        return
    efa = EnhancedFattyAcidMetabolism(smiles)
    c_alpha, c_beta, c_gamma = efa.get_alpha_beta_gamma_carbons(efa.molecule)
    if c_alpha is not None:
        assert isinstance(c_alpha, int)
    if c_beta is not None:
        assert isinstance(c_beta, int)
    if c_gamma is not None:
        assert isinstance(c_gamma, int)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_has_alpha_beta_double_bond(name, smiles):
    if name == "not_fatty_acid":
        return
    efa = EnhancedFattyAcidMetabolism(smiles)
    result = efa.has_alpha_beta_double_bond(efa.molecule)
    assert isinstance(result, bool)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_run_reaction_mechanism(name, smiles):
    if name == "not_fatty_acid":
        return

    efa = EnhancedFattyAcidMetabolism(smiles)

    try:
        product = efa.run_reaction("[CH2:1][CH2:2] >> [CH:1]=[CH:2]", efa.molecule)
        assert product is None or isinstance(product, Chem.Mol)
    except RuntimeError as e:
        # This is expected for some SMILES where the SMARTS doesn't match
        print(f"Expected failure for {name}: {e}")

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_handle_cis_double_bonds(name, smiles):
    if name == "not_fatty_acid":
        return
    efa = EnhancedFattyAcidMetabolism(smiles)
    modified = efa.handle_cis_double_bonds(efa.molecule)
    assert isinstance(modified, Chem.Mol)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_double_bond_info(name, smiles):
    if name == "not_fatty_acid":
        return
    efa = EnhancedFattyAcidMetabolism(smiles)
    info = efa.get_double_bond_info(efa.molecule)
    assert isinstance(info, list)
    for bond in info:
        assert "atoms" in bond
        assert isinstance(bond["atoms"], tuple)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_find_carbonyl_carbon(name, smiles):
    efa = None
    try:
        efa = EnhancedFattyAcidMetabolism(smiles)
    except ValueError:
        return
    idx = efa.find_carbonyl_carbon(efa.molecule)
    if idx is not None:
        assert isinstance(idx, int)

def test_logging_of_reactions():
    efa = EnhancedFattyAcidMetabolism(TEST_MOLECULES["saturated_even"])
    efa.activate_fatty_acid()
    efa.carnitine_shuttle(efa.molecule)
    assert len(efa.reaction_descriptions) >= 1
    assert any("acyl-coa" in desc.lower() for desc in efa.reaction_descriptions)