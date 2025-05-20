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
    "short_chain": "CCCC(=O)O",
}

# SMARTS patterns
CARBOXYLIC_ACID_SMARTS = "[C](=[O])[O]"  # carboxylic acid
THIOESTER_SMARTS = "[C](=[O])[S]"               # thioester

# ------------------ Initialisation, Activation and Carnitine Shuttle ------------------

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_initialization(name, smiles):
    if name == "not_fatty_acid":
        with pytest.raises(ValueError):
            EnhancedFattyAcidMetabolism(smiles)
    else:
        efa = EnhancedFattyAcidMetabolism(smiles)
        assert isinstance(efa.molecule, Chem.Mol)

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_activation_and_transport_detailed(name, smiles):
    if name == "not_fatty_acid":
        return
    
    efa = EnhancedFattyAcidMetabolism(smiles)
    mol = efa.molecule

    # Confirm original molecule has carboxylic acid group (unless test molecule known otherwise)
    ca_pattern = Chem.MolFromSmarts(CARBOXYLIC_ACID_SMARTS)
    
    activated = efa.activate_fatty_acid()
    assert isinstance(activated, Chem.Mol)

    # Activated molecule should NOT have carboxylic acid
    has_ca_after = activated.HasSubstructMatch(ca_pattern)
    assert not has_ca_after, "Activated molecule should no longer have carboxylic acid"

    # Activated molecule SHOULD have thioester group
    te_pattern = Chem.MolFromSmarts(THIOESTER_SMARTS)
    has_te_after = activated.HasSubstructMatch(te_pattern)
    assert has_te_after, "Activated molecule should have thioester group"

    transported = efa.carnitine_shuttle(activated)
    assert isinstance(transported, Chem.Mol)

    # Transported molecule should be identical (no structural change)
    assert transported.GetNumAtoms() == activated.GetNumAtoms()
    assert transported.GetNumBonds() == activated.GetNumBonds()
    assert Chem.MolToSmiles(transported) == Chem.MolToSmiles(activated)

    # Check reaction logging
    assert "Activation" in efa.reaction_steps
    assert "Transport" in efa.reaction_steps

# ------------------ Handling Double Bonds Between Cα-Cβ or Cβ-Cγ ------------------

@pytest.mark.parametrize("name, smiles", TEST_MOLECULES.items())
def test_find_carbonyl_carbon_detailed(name, smiles):
    try:
        efa = EnhancedFattyAcidMetabolism(smiles)
    except ValueError:
        return  # skip invalid molecules

    mol = efa.molecule
    idx = efa.find_carbonyl_carbon(mol)

    # Check if the molecule has a thioester group
    if idx is None:
        return

    # If an index is found, check if it's a valid carbonyl carbon in thioester context
    atom = mol.GetAtomWithIdx(idx)
    assert atom.GetSymbol() == "C", f"Expected carbon at index {idx}"

    # Check neighbors: must include one sulfur (S)
    neighbors = atom.GetNeighbors()
    sulfur_neighbors = [a for a in neighbors if a.GetSymbol() == "S"]
    assert len(sulfur_neighbors) > 0, "Carbonyl carbon must be bonded to sulfur"

    # Check neighbors: must include one oxygen (O) with double bond
    oxygen_neighbors = [a for a in neighbors if a.GetSymbol() == "O"]
    double_bonded_oxygen = False
    for o in oxygen_neighbors:
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx())
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonded_oxygen = True
            break
    assert double_bonded_oxygen, "Carbonyl carbon must have a double bond to oxygen"

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
def test_alpha_beta_gamma_connectivity(name, smiles):

    def is_connected(mol, idx1, idx2):
        return mol.GetBondBetweenAtoms(idx1, idx2) is not None

    if name in ("not_fatty_acid", "no_carbonyl"):
        return
    
    efa = EnhancedFattyAcidMetabolism(smiles)
    mol = efa.molecule
    c_alpha, c_beta, c_gamma = efa.get_alpha_beta_gamma_carbons(mol)

    # If found, alpha must be bonded to carbonyl
    carbonyl_c = efa.find_carbonyl_carbon(mol)
    if c_alpha is not None:
        assert is_connected(mol, c_alpha, carbonyl_c), f"{name}: Cα not connected to carbonyl"
    if c_beta is not None:
        assert is_connected(mol, c_beta, c_alpha), f"{name}: Cβ not connected to Cα"
    if c_gamma is not None:
        assert is_connected(mol, c_gamma, c_beta), f"{name}: Cγ not connected to Cβ"

EXPECTED_ALPHA_BETA_DOUBLE_BOND = {
    "CCCCCC=CC(=O)O": True,
    "CCCCC=CCC(=O)O": False
}

EXPECTED_BETA_GAMMA_DOUBLE_BOND = {
    "CCCCCC=CC(=O)O": False,
    "CCCCC=CCC(=O)O": True
}

@pytest.mark.parametrize("smiles, expected", EXPECTED_ALPHA_BETA_DOUBLE_BOND.items())
def test_has_alpha_beta_double_bond(smiles, expected):
    
    efa = EnhancedFattyAcidMetabolism(smiles)
    activated_mol = efa.activate_fatty_acid()
    
    result = efa.has_alpha_beta_double_bond(activated_mol)

    assert isinstance(result, bool), f"{smiles}: result should be bool"
    assert result == expected, f"{smiles}: expected {expected}, got {result}"

@pytest.mark.parametrize("smiles, expected", EXPECTED_BETA_GAMMA_DOUBLE_BOND.items())
def test_has_beta_gamma_double_bond(smiles, expected):
    
    efa = EnhancedFattyAcidMetabolism(smiles)
    activated_mol = efa.activate_fatty_acid()
    
    result = efa.has_beta_gamma_double_bond(activated_mol)

    assert isinstance(result, bool), f"{smiles}: result should be bool"
    assert result == expected, f"{smiles}: expected {expected}, got {result}"

# --------------------------- run_reaction(self, reaction_smarts, molecule) ---------------------------

SUCCESS_SMILES = "CCCC"          # contains –CH2–CH2–
FAILURE_SMILES = "C=C"           # no –CH2–CH2–
SIMPLE_SMARTS = "[CH2:1][CH2:2] >> [CH:1]=[CH:2]"

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

@pytest.mark.parametrize("smiles, should_match", [
    (SUCCESS_SMILES, True),
    (FAILURE_SMILES, False),
])
def test_run_reaction_matches_or_not(smiles, should_match):
    efa = EnhancedFattyAcidMetabolism(smiles)
    product = efa.run_reaction(SIMPLE_SMARTS, efa.molecule)

    if should_match:
        # We expect a product and a new double bond
        assert isinstance(product, Chem.Mol)
        assert "=" in Chem.MolToSmiles(product), "double bond not created"
    else:
        # We expect no product
        assert product is None

def test_run_reaction_invalid_smarts():
    efa = EnhancedFattyAcidMetabolism("CC")
    # "Failed to run reaction" string appears in the error message raised by "run_reaction"
    with pytest.raises(RuntimeError, match="Failed to run reaction"):
        efa.run_reaction("this_is_not_smarts", efa.molecule)

# --------------------------- beta_oxidation_cycle(self, acyl_coa_mol) ---------------------------

def test_beta_oxidation_cycle_basic():
    smiles = "CCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    acyl_coa = efa.activate_fatty_acid()
    
    product = efa.beta_oxidation_cycle(acyl_coa)

    assert isinstance(product, Chem.Mol)
    assert efa.atp_yield > 0  # Should gain ATP
    assert "Dehydrogenation" in efa.reaction_steps or "Isomerase‑assisted dehydrogenation" in efa.reaction_steps
    assert len(efa.cycle_log) == 1
    assert "steps" in efa.cycle_log[0]
    assert len(efa.cycle_log[0]["steps"]) >= 1

def test_beta_oxidation_cycle_with_beta_gamma_double_bond():
    smiles = "CC=CCC(=O)O"  # β-γ double bond should trigger isomerase-assisted path
    efa = EnhancedFattyAcidMetabolism(smiles)
    acyl_coa = efa.activate_fatty_acid()

    product = efa.beta_oxidation_cycle(acyl_coa)

    assert isinstance(product, Chem.Mol)
    assert "Isomerase‑assisted dehydrogenation" in efa.reaction_steps
    assert any("Isomerase" in desc for desc in efa.reaction_descriptions)

def test_beta_oxidation_cycle_logs_correctly():
    smiles = "CCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    acyl_coa = efa.activate_fatty_acid()

    efa.beta_oxidation_cycle(acyl_coa)

    # Check ATP yield
    assert efa.atp_yield in (4.0, 1.5 + 2.5)  # FADH2 + NADH

    # Check reaction steps and logs
    assert len(efa.reaction_steps) >= 2
    assert all(isinstance(step, str) for step in efa.reaction_steps)
    assert all(isinstance(res, Chem.Mol) for res in efa.reaction_results)
    assert all("smarts" in s for s in efa.cycle_log[0]["steps"])

# ------------------------------ run_complete_oxidation(self) ------------------------------

def test_run_complete_oxidation_even_chain():
    # Example: octanoic acid (8 carbons, even chain)
    smiles = "CCCCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    
    result = efa.run_complete_oxidation()
    
    assert isinstance(result, dict)
    assert result['final_products']['acetyl_coa_count'] > 0
    assert not result['final_products']['propionyl_coa']  # even chain → no propionyl-CoA
    assert isinstance(result['reaction_steps'], list)
    assert isinstance(result['reaction_results'], list)
    assert isinstance(result['reaction_descriptions'], list)
    assert result['total_atp_yield'] >= 0

def test_run_complete_oxidation_odd_chain():
    # Example: nonanoic acid (9 carbons, odd chain)
    smiles = "CCCCCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    
    result = efa.run_complete_oxidation()
    
    assert isinstance(result, dict)
    assert result['final_products']['acetyl_coa_count'] > 0
    assert result['final_products']['propionyl_coa']  # odd chain → propionyl-CoA present
    assert isinstance(result['reaction_steps'], list)
    assert isinstance(result['reaction_results'], list)
    assert isinstance(result['reaction_descriptions'], list)
    assert result['total_atp_yield'] >= 0

def test_run_complete_oxidation_short_chain_edge_case():
    # Example: acetic acid (2 carbons, chain too short to beta-oxidize)
    smiles = "CC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    
    result = efa.run_complete_oxidation()
    
    # Should handle edge case gracefully
    assert isinstance(result, dict)
    assert result['final_products']['acetyl_coa_count'] == 1 or result['final_products']['propionyl_coa'] is False
    assert 'Error: Chain too short to process' in efa.reaction_cycles or len(efa.reaction_cycles) > 0
    assert result['total_atp_yield'] >= 0

def test_run_complete_oxidation_activation_failure(monkeypatch):
    # Simulate failure in activation step by monkeypatching activate_fatty_acid
    smiles = "CCCCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    
    # Monkeypatch "activate_fatty_acid" to raise error on call
    # Replaces the activate_fatty_acid method on the efa instance with fail_activation, which always raises an error.
    # When "run_complete_oxidation()" is called, instead of running the real activation method, it runs the fake one that triggers a failure.
    # Tests how "run_complete_oxidation" behaves when activation fails

    def fail_activation():
        raise RuntimeError("Simulated activation failure")

    monkeypatch.setattr(efa, 'activate_fatty_acid', fail_activation)
    
    with pytest.raises(RuntimeError, match="Simulated activation failure"):
        efa.run_complete_oxidation()

def test_run_complete_oxidation_beta_oxidation_cycle_failure(monkeypatch):
    # Simulate failure in beta_oxidation_cycle to test error propagation
    smiles = "CCCCCCCC(=O)O"
    efa = EnhancedFattyAcidMetabolism(smiles)
    
    # Allow activation and transport to succeed normally
    # Monkeypatch beta_oxidation_cycle to raise error on call

    def fail_beta_oxidation_cycle(mol):
        raise RuntimeError("Simulated beta-oxidation cycle failure")
    
    monkeypatch.setattr(efa, 'beta_oxidation_cycle', fail_beta_oxidation_cycle)
    
    with pytest.raises(RuntimeError, match="Simulated beta-oxidation cycle failure"):
        efa.run_complete_oxidation()

# --------------------------------- prepare_data_for_visualization(self) ---------------------------------

def test_prepare_data_for_visualization_basic():
    smiles = "CCO"
    efa = EnhancedFattyAcidMetabolism(smiles) 
    
    # Setup example data
    mol = Chem.MolFromSmiles("CCO")  # Ethanol molecule
    efa.reaction_steps = ["Step 1"]
    efa.reaction_descriptions = ["Description of step 1"]
    efa.reaction_results = [mol]
    efa.atp_yield = 10.0
    efa.cycle_log = [{"cycle_number": 1, "steps": []}]

    result = efa.prepare_data_for_visualization()

    assert result["metadata"]["total_steps"] == 1
    assert result["metadata"]["total_atp_yield"] == 10.0
    assert result["metadata"]["total_cycles"] == 1

    step = result["steps"][0]
    assert step["index"] == 0
    assert step["step_number"] == 1
    assert step["name"] == "Step 1"
    assert step["description"] == "Description of step 1"
    assert step["smiles"] == Chem.MolToSmiles(mol, canonical=True)
    assert step["formula"] == Chem.rdMolDescriptors.CalcMolFormula(mol)

    assert result["cycles"] == efa.cycle_log


def test_prepare_data_for_visualization_invalid_molecule(monkeypatch):
    smiles = "CCO"
    efa = EnhancedFattyAcidMetabolism(smiles) 

    # Create a dummy molecule that fails sanitization
    mol = Chem.MolFromSmiles("CCO")
    efa.reaction_steps = ["Step 1"]
    efa.reaction_descriptions = ["Description of step 1"]
    efa.reaction_results = [mol]
    efa.atp_yield = 5.0
    efa.cycle_log = []

    # Monkeypatch Chem.SanitizeMol to raise an exception to simulate invalid molecule
    def raise_error(m):
        raise ValueError("Invalid molecule")
    monkeypatch.setattr(Chem, "SanitizeMol", raise_error)

    result = efa.prepare_data_for_visualization()

    step = result["steps"][0]
    assert step["smiles"] == "Invalid"
    assert step["formula"] == "Invalid"

def test_prepare_data_for_visualization_empty_data():
    smiles = "CCO"
    efa = EnhancedFattyAcidMetabolism(smiles)

    efa.reaction_steps = []
    efa.reaction_descriptions = []
    efa.reaction_results = []
    efa.atp_yield = 0
    efa.cycle_log = []

    result = efa.prepare_data_for_visualization()

    assert result["metadata"]["total_steps"] == 0
    assert result["metadata"]["total_atp_yield"] == 0
    assert result["metadata"]["total_cycles"] == 0
    assert result["steps"] == []
    assert result["cycles"] == []