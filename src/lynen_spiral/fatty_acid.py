"""Module for fatty acid metabolism analysis."""
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import IPythonConsole
import matplotlib.pyplot as plt
from enum import Enum
from typing import Dict, List, Tuple

class FattyAcidType(Enum):
    SATURATED = "saturated"
    MONOUNSATURATED = "monounsaturated"
    POLYUNSATURATED = "polyunsaturated"
    
class FattyAcidMetabolism:
    def __init__(self, smiles: str):
        """
        Initialize fatty acid metabolism analyzer.
        
        Parameters
        ----------
        smiles : str
            SMILES representation of the fatty acid
        """
        self.molecule = Chem.MolFromSmiles(smiles)
        if not self.molecule:
            raise ValueError("Invalid SMILES string")
            
        # Basic properties
        self.chain_length = self._get_carbon_chain_length()
        self.double_bonds = self._count_double_bonds()
        self.double_bond_positions = self._get_double_bond_positions()
        self.type = self._determine_fatty_acid_type()
        
    def _get_carbon_chain_length(self) -> int:
        """Calculate the length of the carbon chain."""
        return len([atom for atom in self.molecule.GetAtoms() 
                   if atom.GetSymbol() == 'C'])
    
    def _count_double_bonds(self) -> int:
        """Count the number of C=C double bonds."""
        pattern = Chem.MolFromSmarts('C=C')
        return len(self.molecule.GetSubstructMatches(pattern))
    
    def _get_double_bond_positions(self) -> List[int]:
        """Get the positions of double bonds."""
        pattern = Chem.MolFromSmarts('C=C')
        matches = self.molecule.GetSubstructMatches(pattern)
        return [match[0] for match in matches]
    
    def _determine_fatty_acid_type(self) -> FattyAcidType:
        """Determine the type of fatty acid."""
        if self.double_bonds == 0:
            return FattyAcidType.SATURATED
        elif self.double_bonds == 1:
            return FattyAcidType.MONOUNSATURATED
        else:
            return FattyAcidType.POLYUNSATURATED
    
    def calculate_activation_energy(self) -> Dict[str, float]:
        """
        Calculate ATP cost for fatty acid activation.
        
        The activation of fatty acids requires 2 ATP:
        1. ATP → AMP + PPi (equivalent to 2 ATP)
        """
        return {
            'ATP_cost': 2,
            'description': 'Conversion to acyl-CoA via acyl-CoA synthetase'
        }
    
    def calculate_atp_yield(self) -> dict:
        """
        Calculate ATP yield from complete oxidation.
        
        Returns
        -------
        dict
            Breakdown of ATP production steps
        """
        # First calculate activation cost
        activation_cost = self.calculate_activation_energy()['ATP_cost']
        
        # Number of β-oxidation cycles
        cycles = (self.chain_length - 2) // 2
        
        # Base yields per cycle
        fadh2_per_cycle = 1  # FADH2 → 1.5 ATP
        nadh_per_cycle = 1   # NADH → 2.5 ATP
        acetyl_coa_per_cycle = 1  # Each acetyl-CoA → 10 ATP via citric acid cycle
        
        # Additional NADH from unsaturated FA metabolism
        extra_nadh = 0
        if self.type != FattyAcidType.SATURATED:
            # Each double bond requires one additional NADH
            extra_nadh = self.double_bonds
        
        # Handle odd-chain fatty acids
        is_odd_chain = self.chain_length % 2 != 0
        propionyl_coa_atp = 0
        if is_odd_chain:
            # Propionyl-CoA metabolism yields additional ATP
            propionyl_coa_atp = 13  # Simplified value
        
        # Calculate totals
        total_fadh2 = cycles * fadh2_per_cycle
        total_nadh = (cycles * nadh_per_cycle) + extra_nadh
        total_acetyl_coa = cycles * acetyl_coa_per_cycle
        
        if is_odd_chain:
            total_acetyl_coa -= 1  # Last fragment is propionyl-CoA instead
        
        # Calculate final ATP yield
        total_atp = (
            (total_fadh2 * 1.5) +  # FADH2 contribution
            (total_nadh * 2.5) +   # NADH contribution
            (total_acetyl_coa * 10) +  # Acetyl-CoA contribution
            propionyl_coa_atp -    # Propionyl-CoA if odd-chain
            activation_cost        # Subtract activation cost
        )
        
        return {
            'activation_cost': activation_cost,
            'β-oxidation_cycles': cycles,
            'FADH2_ATP': total_fadh2 * 1.5,
            'NADH_ATP': total_nadh * 2.5,
            'acetyl_CoA_ATP': total_acetyl_coa * 10,
            'propionyl_CoA_ATP': propionyl_coa_atp,
            'is_odd_chain': is_odd_chain,
            'unsaturation_type': self.type.value,
            'double_bonds': self.double_bonds,
            'total_ATP': total_atp
        }
    
    def visualize_steps(self):
        """Visualize the β-oxidation steps including activation."""
        steps = []
        descriptions = []
        
        # Add initial molecule
        steps.append(self.molecule)
        descriptions.append("Initial fatty acid")
        
        # Add activation step (conversion to acyl-CoA)
        activated_smiles = self._simulate_activation(Chem.MolToSmiles(self.molecule))
        activated_mol = Chem.MolFromSmiles(activated_smiles)
        steps.append(activated_mol)
        descriptions.append("Activated fatty acid (acyl-CoA)")
        
        current_mol = activated_mol
        cycle = 1
        
        while len([atom for atom in current_mol.GetAtoms() 
                  if atom.GetSymbol() == 'C']) > 2:
            
            # Handle unsaturated bonds if present
            if self.type != FattyAcidType.SATURATED and cycle in self._get_double_bond_positions():
                # Add isomerization step
                isomerized_smiles = self._simulate_isomerization(Chem.MolToSmiles(current_mol))
                isomerized_mol = Chem.MolFromSmiles(isomerized_smiles)
                steps.append(isomerized_mol)
                descriptions.append(f"Cycle {cycle}: Isomerization")
                current_mol = isomerized_mol
            
            # Regular β-oxidation steps
            oxidized_smiles = self._remove_two_carbons(Chem.MolToSmiles(current_mol))
            current_mol = Chem.MolFromSmiles(oxidized_smiles)
            steps.append(current_mol)
            descriptions.append(f"Cycle {cycle}: β-oxidation")
            cycle += 1
        
        # Visualization
        fig, axes = plt.subplots(len(steps), 1, figsize=(10, 5*len(steps)))
        if len(steps) == 1:
            axes = [axes]
            
        for idx, (mol, desc) in enumerate(zip(steps, descriptions)):
            img = Draw.MolToImage(mol)
            axes[idx].imshow(img)
            axes[idx].axis('off')
            axes[idx].set_title(f'Step {idx+1}: {desc}')
        
        plt.tight_layout()
        return fig
    
    def _simulate_activation(self, smiles: str) -> str:
        """Simulate conversion to acyl-CoA."""
        # This is a simplified representation - in reality would be more complex
        return smiles.replace("O", "SCoA")
    
    def _simulate_isomerization(self, smiles: str) -> str:
        """Simulate isomerization of double bonds."""
        # Simplified - would need more complex logic for real implementation
        return smiles
    
    def _remove_two_carbons(self, smiles: str) -> str:
        """Simulate removal of two carbons."""
        mol = Chem.MolFromSmiles(smiles)
        atoms = mol.GetAtoms()
        n_carbons = len([atom for atom in atoms if atom.GetSymbol() == 'C'])
        return f"C{'C' * (n_carbons-4)}(=O)SCoA"  # Simplified
