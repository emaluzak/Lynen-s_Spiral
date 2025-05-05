import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, MolFromSmiles, Draw

import sys
import os
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
sys.path.append(os.path.abspath(r"C:\Users\ipeki\git\Lynen-s_Spiral\Lynen-s_Spiral\src\lynen_spiral"))

from enhanced_fatty_acid import EnhancedFattyAcidMetabolism #type: ignore

# Test cases for input handling
@pytest.mark.parametrize("input_value, expected_smiles", [
    ("CCCCCCCCCCCCCCCC(=O)O", "CCCCCCCCCCCCCCCC(=O)O"),  # SMILES input
    ("palmitic acid", "CCCCCCCCCCCCCCCC(=O)O"),          # Common name input
    ("C16:0", "CCCCCCCCCCCCCCCC(=O)O"),                  # CX:Y notation (saturated)
    ("C18:1(9)", "CCCCCCCCC=CCCCCCCC(=O)O"),             # CX:Y notation (monounsaturated)
    ("C18:2(9,12)", "CCCCCC=CCC=CCCCCCC(=O)O"),          # CX:Y notation (polyunsaturated)
])
def test_input_handling(input_value, expected_smiles):
    efa = EnhancedFattyAcidMetabolism(input_value)
    parsed_smiles = Chem.MolToSmiles(efa.molecule)
    print(f"Input: {input_value}")
    print(f"Expected: {expected_smiles}")
    print(f"Generated: {parsed_smiles}")
    assert parsed_smiles == expected_smiles

# Test cases for reaction SMARTS 1
@pytest.mark.parametrize("reaction_name, smarts, input_smiles, expected_product", [
    ("dehydrogenation", "[CH2:1][CH2:2][C:3](=O)[S:4] >> [CH:1]=[CH:2][C:3](=O)[S:4]", 
     "CCCC(=O)S", "CC=CC(=O)S"), # Dehydrogenation

    ("hydration", "[CH:1]=[CH:2][C:3](=O)[S:4] >> [CH:1]([OH])[CH2:2][C:3](=O)[S:4]", 
     "CC=CC(=O)S", "CC(O)CC(=O)S"),  # Hydration

    ("oxidation", "[CH:1]([OH])[CH2:2][C:3](=O)[S:4] >> [C:1](=O)[CH2:2][C:3](=O)[S:4]", 
     "CC(O)CC(=O)S", "CC(=O)CC(=O)S"), # Oxidation

    ("thiolysis", "[CH3:1][C:2](=O)[CH2:3][C:4](=O)[S:5] >> [CH3:1][C:2](=O)[S:5]", 
     "CC(=O)CC(=O)S", "CC(=O)S") # Thiolysis
])
def test_reaction_smarts(reaction_name, smarts, input_smiles, expected_product):
    efa = EnhancedFattyAcidMetabolism(input_smiles)
    reactant = efa.molecule
    product = efa.run_reaction(smarts, reactant)
    
    # Make sure to access the first product if multiple are returned
    if product:
        product_mol = product[0] if isinstance(product, tuple) else product
         # Remove explicit hydrogens
        product_mol = Chem.RemoveHs(product_mol)
        expected_mol = Chem.MolFromSmiles(expected_product)
        expected_mol = Chem.RemoveHs(expected_mol)
        assert Chem.MolToSmiles(product_mol, canonical=True, isomericSmiles=False) == Chem.MolToSmiles(MolFromSmiles(expected_product), canonical=True, isomericSmiles=False)
    else:
        assert False, "No product generated"


# Test cases for reaction SMARTS 2
@pytest.mark.parametrize("reaction_name, smarts, input_smiles, expected_product", [
    # Secondary alcohol oxidation to ketone (with CoA activation)
    ("oxidation", "[CH:1]([OH])[CH2:2][C:3](=O)[S:4] >> [C:1](=O)[CH2:2][C:3](=O)[S:4]",
     "CC(O)CC(=O)S", "CC(=O)CC(=O)S"), # Oxidation

    ("hydration", "[CH:1]=[CH:2][C:3](=O)[S:4] >> [CH:1]([OH])[CH2:2][C:3](=O)[S:4]", 
     "CC=CC(=O)S", "CC(O)CC(=O)S"), # Hydration

    # Dehydrogenation (double bond formation near CoA group)
    ("dehydrogenation", "[CH2:1][CH2:2][C:3](=O)[S:4] >> [CH:1]=[CH:2][C:3](=O)[S:4]", 
     "CCCC(=O)S", "CC=CC(=O)S"), # Dehydrogenation

    # Thiolysis (cleaves β-keto thioester into a shorter thioester + acid)
    ("thiolysis", "[CH3:1][C:2](=O)[CH2:3][C:4](=O)[S:5] >> [CH3:1][C:2](=O)[S:5]",
     "CC(=O)CC(=O)S", "CC(=O)S"), # Thiolysis
])
def test_fatty_acid_reactions(reaction_name, smarts, input_smiles, expected_product):
    efa = EnhancedFattyAcidMetabolism(input_smiles)
    reactant = efa.molecule
    product = efa.run_reaction(smarts, reactant)

    if product:
        product_mol = product[0] if isinstance(product, tuple) else product
        product_mol = Chem.RemoveHs(product_mol)
        expected_mol = Chem.MolFromSmiles(expected_product)
        expected_mol = Chem.RemoveHs(expected_mol)
        assert Chem.MolToInchiKey(product_mol) == Chem.MolToInchiKey(expected_mol)

    else:
        assert False, "No product generated"

# Test cases for complete oxidation (beta-oxidation cycles)
@pytest.mark.parametrize("input_value, expected_atp_yield", [
    ("CCCCCCCCCCCCCCCC(=O)O", 129),  # Palmitic acid (C16:0)
    ("CCCCCCCCCCCCCCCCCC(=O)O", 163),  # Stearic acid (C18:0)
    ("CCCCCCCCC=CCCCCCCC(=O)O", 120),  # Oleic acid (C18:1)
    ("CCCCCC=CCC=CCCCCCC(=O)O", 114),  # Linoleic acid (C18:2)
    ("CCC=CCC=CCC=CCCCCCC(=O)O", 108),  # Alpha-linolenic acid (C18:3)
])
def test_complete_oxidation(input_value, expected_atp_yield):
    efa = EnhancedFattyAcidMetabolism(input_value)
    results = efa.run_complete_oxidation()
    assert results['total_atp_yield'] == expected_atp_yield


# Test error handling for invalid inputs
@pytest.mark.parametrize("invalid_input", [
    "InvalidSMILES",
    "C18:2(9,12,15,18)",  # Too many double bonds
    "Unknown fatty acid",
])
def test_invalid_inputs(invalid_input):
    with pytest.raises(ValueError):
        EnhancedFattyAcidMetabolism(invalid_input)


# Test error handling for failed reactions
def test_reaction_produces_no_products():
    efa = EnhancedFattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")  # Palmitic acid
    invalid_smarts = "[C:1 >> [C:1]="  # malformed — will fail to compile
    with pytest.raises(RuntimeError):
        efa.run_reaction(invalid_smarts, efa.molecule)


# Test data structure for visualization
def test_prepare_data_for_visualization():
    efa = EnhancedFattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")  # Palmitic acid
    results = efa.run_complete_oxidation()
    data = efa.prepare_data_for_visualization()
    
    # Check overall structure
    assert "steps" in data
    assert "total_atp_yield" in data
    assert data["total_atp_yield"] > 0

    # Check each step
    for step in data["steps"]:
        assert "step_number" in step
        assert "step_name" in step
        assert "description" in step
        assert "molecule_smiles" in step
        assert "formula" in step


def test_prepare_data_for_visualization_with_image():
    efa = EnhancedFattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")  # Palmitic acid
    results = efa.run_complete_oxidation()
    data = efa.prepare_data_for_visualization()
    
    # Check overall structure
    assert "steps" in data
    assert "total_atp_yield" in data
    assert data["total_atp_yield"] > 0

    # Check each step and validate molecules
    for step in data["steps"]:
        assert "step_number" in step
        assert "step_name" in step
        assert "description" in step
        assert "molecule_smiles" in step
        assert "formula" in step

    # Now let's visualize just the first molecule
    first_step = data["steps"][12]  # Get the twelvth step
    molecule_smiles = first_step["molecule_smiles"]  # Get the SMILES of the first molecule
    molecule = Chem.MolFromSmiles(molecule_smiles)  # Convert SMILES to RDKit molecule object

    if molecule is not None:
        # Visualize the first molecule
        img = Draw.MolToImage(molecule)
        img.show()  # Show the image of the first molecule
    else:
        print(f"Error: Invalid molecule at step 1 with SMILES: {molecule_smiles}")


# Test the visualize_reaction_sequence method
def test_visualize_reaction_sequence():
    efa = EnhancedFattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")  # Palmitic acid
    results = efa.run_complete_oxidation()  # Run the oxidation steps
    
    # Now test if visualization is working properly
    fig = efa.visualize_reaction_sequence()  # Assuming it returns a matplotlib figure
    
    # Check if the figure was created (i.e., not None)
    assert fig is not None, "Visualization figure is None"
    
    # Optionally, you can check if the number of steps is correct
    # Let's assume you expect 5 steps (or adjust based on your real output)
    assert len(fig.axes) == len(efa.reaction_steps), "Number of steps in visualization does not match reaction steps"
    
    # Check if each step has a title (assuming you set titles in the visualization)
    for ax in fig.axes:
        assert ax.get_title() != "", "Step title is empty"

    # Check if the plot is created (in some environments you might want to check if it can be saved or shown)
    # You can also assert based on the number of images or any other specific visualization-related assertions
    
    # Optionally, save the figure to test saving
    fig.savefig(r"C:\Users\ipeki\OneDrive\Masaüstü\lynen test saves\test_visualization_output.png")

    # Check if the saved figure exists (optional)
    import os
    assert os.path.exists("test_visualization_output.png"), "Visualization output file was not saved"

def test_prepare_data_for_visualization(): # this test is probably wrong
    # Arrange
    mol = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)O")  # Palmitic acid
    pathway = EnhancedFattyAcidMetabolism()
    pathway.reaction_steps = ["Fatty Acid Activation"]
    pathway.reaction_descriptions = ["Activation of palmitic acid to acyl-CoA"]
    pathway.reaction_results = [mol]
    pathway.atp_yield = 129

    # Act
    data = pathway.prepare_data_for_visualization()

    # Assert
    assert data["metadata"]["total_steps"] == 1
    assert data["metadata"]["total_atp_yield"] == 129

    step = data["steps"][0]
    assert step["name"] == "Fatty Acid Activation"
    assert step["smiles"] == "CCCCCCCCCCCCCCCC(=O)O"
    assert step["formula"] == "C16H32O2"
    assert step["step_number"] == 1
