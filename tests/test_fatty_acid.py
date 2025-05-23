import pytest
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import sys
import os
sys.path.append(os.path.abspath(r"src\lynen_spiral\Lynen_spiral_visualisation"))
from fatty_acid import FattyAcidType, FattyAcidMetabolism

def test_saturated_even_chain():
    # Test with palmitic acid (C16:0)
    palmitic = FattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")
    atp_yield = palmitic.calculate_atp_yield()
    
    assert atp_yield['β-oxidation_cycles'] == 7
    assert atp_yield['is_odd_chain'] is False
    assert atp_yield['total_ATP'] == 106

def test_unsaturated_fatty_acid():
    # Test with oleic acid (C18:1)
    oleic = FattyAcidMetabolism("CCCCCCCC=CCCCCCCCC(=O)O")
    atp_yield = oleic.calculate_atp_yield()
    
    assert atp_yield['double_bonds'] == 1
    assert atp_yield['unsaturation_type'] == FattyAcidType.MONOUNSATURATED.value

def test_odd_chain_fatty_acid():
    # Test with C15:0
    pentadecanoic = FattyAcidMetabolism("CCCCCCCCCCCCCCC(=O)O")
    atp_yield = pentadecanoic.calculate_atp_yield()
    
    assert atp_yield['is_odd_chain'] is True
    assert atp_yield['propionyl_CoA_ATP'] > 0

def test_invalid_smiles():
    with pytest.raises(ValueError):
        FattyAcidMetabolism("invalid_smiles")