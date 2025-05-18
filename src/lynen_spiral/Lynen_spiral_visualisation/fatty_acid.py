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
        
        Parameters:
            smiles : str
                     SMILES representation of the fatty acid
        """

        self.molecule = Chem.MolFromSmiles(smiles)
        if not self.molecule:
            raise ValueError("Invalid SMILES string")
            
        # Basic properties
        self.chain_length = self._get_carbon_chain_length(self.molecule)
        self.double_bonds = self._count_double_bonds()
        self.double_bond_positions = self._get_double_bond_positions()
        self.type = self._determine_fatty_acid_type()
    
    def _get_carbon_chain_length(self, mol) -> int:

        """
        Calculate the total number of carbon atoms in the given molecule.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object.

        Returns:
            int: Total number of carbon atoms in the molecule.
        """
        
        return len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'])

    
    def _count_double_bonds(self) -> int:

        """
        Count the number of carbon–carbon double bonds (C=C) in the molecule 
        using a SMARTS pattern.

        Returns:
            int: The number of C=C double bonds found in the molecule.
        """

        pattern = Chem.MolFromSmarts('C=C')
        return len(self.molecule.GetSubstructMatches(pattern))
    
    def _get_double_bond_positions(self) -> List[int]:

        """
        Identify the positions of double bonds in the fatty acid 
        molecule using a SMART pattern.

        Returns:
            List[int]: A list of atom indices corresponding to the first atom 
            in each C=C double bond.
        """

        pattern = Chem.MolFromSmarts('C=C')
        matches = self.molecule.GetSubstructMatches(pattern)
        return [match[0] for match in matches]
    
    def _determine_fatty_acid_type(self) -> FattyAcidType:

        """
        Classify the fatty acid based on the number of double bonds.

        Returns:
            FattyAcidType: One of the following enum values indicating the saturation level:
                - SATURATED: No double bonds
                - MONOUNSATURATED: Exactly one double bond
                - POLYUNSATURATED: More than one double bond
        """

        if self.double_bonds == 0:
            return FattyAcidType.SATURATED
        elif self.double_bonds == 1:
            return FattyAcidType.MONOUNSATURATED
        else:
            return FattyAcidType.POLYUNSATURATED
    
    def calculate_activation_energy(self) -> Dict[str, float]:

        """
        Calculate the ATP cost associated with fatty acid activation.

        Fatty acids are activated to acyl-CoA by the enzyme acyl-CoA synthetase,
        which hydrolyzes ATP to AMP and pyrophosphate (PPi). This reaction is 
        energetically equivalent to consuming 2 ATP molecules.

        Returns:
            Dict[str, float]: A dictionary containing:
                - 'ATP_cost': (float) The energetic cost in ATP equivalents (2 ATP)
                - 'description': (str) A short explanation of the activation step
        """
         
        return {
            'ATP_cost': 2,
            'description': 'Conversion to acyl-CoA via acyl-CoA synthetase'
        }
    
    def calculate_atp_yield(self) -> dict:

        """
        Calculate the total ATP yield from complete β-oxidation of a fatty acid.

        This method computes the ATP contribution from:
        - FADH2 and NADH generated during β-oxidation cycles,
        - Acetyl-CoA entering the citric acid cycle,
        - Adjustments for unsaturated fatty acids (e.g., isomerase/reductase effects),
        - Odd-chain fatty acids generating propionyl-CoA,
        - Subtracting activation energy costs.

        Returns:
            dict: A breakdown of ATP yield components including:
                - Total FADH2 and NADH ATP equivalents
                - Total Acetyl-CoA yield and contribution
                - Extra NADH from unsaturation handling
                - ATP penalty from NADPH-consuming reductase (if applicable)
                - ATP yield from propionyl-CoA (odd-chain FAs)
                - Activation energy cost
                - Final total ATP yield

        Notes:
            - Assumes 1 FADH2 = 1.5 ATP, 1 NADH = 2.5 ATP, 1 Acetyl-CoA = 10 ATP.
            - Simplifies propionyl-CoA contribution as a fixed ATP gain for odd-chain FAs.
            - May include a reductase penalty for even-position double bonds in unsaturated FAs.
            - Uses "calculate_activation_energy" method.
        """

        # First calculate activation cost
        activation_cost = self.calculate_activation_energy()['ATP_cost']

        # Number of β-oxidation cycles
        cycles = (self.chain_length - 2) // 2

        # Base yields per cycle
        fadh2_per_cycle = 1  # FADH2 → 1.5 ATP
        nadh_per_cycle = 1   # NADH → 2.5 ATP
        acetyl_coa_per_cycle = 1  # Each acetyl-CoA → 10 ATP via citric acid cycle

        # Initialize propionyl_coa_atp here before modifying it in the loop
        propionyl_coa_atp = 0  # This will be modified for odd-chain fatty acids later
    
        # Additional NADH from unsaturated FA metabolism
        extra_nadh = 0
        isomerization_step = None  # Initialize the variable for tracking isomerization step
        if self.type != FattyAcidType.SATURATED:
            for position in self.double_bond_positions:
                if position % 2 == 1:
                    # Enoyl-CoA Isomerase - no extra cost
                    extra_nadh += 1
                    isomerization_step = "Isomerization step completed"
                else:
                    # 2,4-Dienoyl-CoA Reductase - NADPH cost reduces yield
                    propionyl_coa_atp -= 2 / 2.5  # Adjust ATP for reductase step (Convert NADPH cost to ATP units)
                    isomerization_step = "Isomerization step with reductase"

        # Handle odd-chain fatty acids
        is_odd_chain = self.chain_length % 2 != 0
        propionyl_coa_atp = 0
        if is_odd_chain:
            # Propionyl-CoA metabolism yields additional ATP
            propionyl_coa_atp = 15  # Simplified value

        # Calculate totals
        total_fadh2 = cycles * fadh2_per_cycle
        total_nadh = (cycles * nadh_per_cycle) + extra_nadh
        total_acetyl_coa = cycles * acetyl_coa_per_cycle + 1

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
            'total_ATP': total_atp,
            'isomerization_step': isomerization_step  # Add this line to track isomerization
        }
    
    def visualize_steps(self):

        """
        Generate a step-by-step visual representation of the β-oxidation process.

        Simulate and visualize the key stages of fatty acid metabolism:
        1. The initial unmodified fatty acid.
        2. Activation to the acyl-CoA form.
        3. Repeated β-oxidation cycles, including any necessary isomerizations 
           (e.g., for cis double bonds) and two-carbon removals.

        Each step is recorded and rendered as a molecular image using RDKit and 
        displayed with its corresponding description in a vertically stacked matplotlib plot.

        Returns:
            matplotlib.figure.Figure: A matplotlib figure containing molecular images 
            and stepwise annotations of the β-oxidation pathway.
    
        Notes:
            - This visualization is simplified and symbolic. Chemical accuracy is 
              approximated to provide pedagogical clarity.
            - Isomerization steps are simulated only at cycles with known double bond positions.
        """

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
        
        while len([atom for atom in current_mol.GetAtoms() if atom.GetSymbol() == 'C']) > 2:
            # Apply isomerization for double bonds
            if cycle in self.double_bond_positions:
                current_mol = Chem.MolFromSmiles(
                self._simulate_isomerization(Chem.MolToSmiles(current_mol), cycle)
                )
                steps.append(current_mol)
                descriptions.append(f"Cycle {cycle}: Isomerization")

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

        """
        Simulate conversion to acyl-CoA.
        
        Args:
            smiles (str): The SMILES representation of the unactivated fatty acid.

        Returns:
            str: A modified SMILES string representing the acyl-CoA form.

        Note:
            This function does not generate a chemically valid acyl-CoA structure.
            It serves as a symbolic placeholder for illustrative purposes.
        """
        
        # This is a simplified representation - in reality would be more complex
        return smiles.replace("O", "SCoA")
    
    # (Ipek) Modified the _simulate_isomerization method to handle fatty acids of different lengths and return their SMILES.
    def _simulate_isomerization(self, smiles: str, position: int) -> str:
        """
        Simulate enzyme-specific isomerization of double bonds.

        Parameters
        ----------
        smiles : str
            SMILES representation of the molecule.
        position : int
            Position of the double bond to be isomerized.

        Returns
        -------
        str
            Updated SMILES string after isomerization.
        """
        if position % 2 == 1:
            # Odd-numbered double bond - enoyl-CoA isomerase
            print(f"Applying enoyl-CoA isomerase at position {position}")
        else:
            # Even-numbered double bond - 2,4-dienoyl-CoA reductase
            print(f"Applying 2,4-dienoyl-CoA reductase at position {position}")

        # Simplified logic - real implementation would modify the molecule structure.
        return smiles  # Placeholder for actual modification logic

    
    def _remove_two_carbons(self, smiles: str) -> str:
        """
        Simulate the shortening of a fatty acid chain by removing two carbon atoms.

        This method interprets a SMILES string, counts the carbon atoms, and constructs 
        a new simplified thioester SMILES string with two fewer carbons (excluding the 
        carbonyl and CoA moiety which are retained as C(=O)SCoA).

        Parameters:
            smiles (str): The SMILES representation of the input fatty acid molecule.

        Returns:
            str: A simplified SMILES string representing the shortened molecule.
        """
        mol = Chem.MolFromSmiles(smiles)
        atoms = mol.GetAtoms()
        n_carbons = len([atom for atom in atoms if atom.GetSymbol() == 'C'])
        return f"C{'C' * (n_carbons-4)}(=O)SCoA"  # Simplified
    

    # (Ipek) New addition: glycerol
    def metabolize_glycerol(self) -> Dict[str, float]:
      
        """
        Simulate glycerol metabolism.

        Glycerol is converted to DHAP via two enzymatic steps:
        1. Glycerol → Glycerol-3-Phosphate (requires 1 ATP)
            (Glycerol is phosphorylated to glycerol-3-phosphate by the enzyme glycerol kinase.)
        2. Glycerol-3-Phosphate → DHAP (produces 1 NADH)
            (Glycerol-3-phosphate is then oxidized to dihydroxyacetone phosphate (DHAP) by glycerol-3-phosphate dehydrogenase.)
            (DHAP enters glycolysis or gluconeogenesis.)

        Returns
        -------
        dict
            ATP and intermediate yields from glycerol metabolism.
        """

        # Step 1: ATP cost for phosphorylation
        atp_cost = 1

        # Step 2: NADH production from oxidation
        nadh_yield = 1
        atp_from_nadh = nadh_yield * 2.5  # NADH yields ~2.5 ATP

        # Net ATP yield from glycerol metabolism
        net_atp = atp_from_nadh - atp_cost

        return {
            'ATP_cost': atp_cost,
            'NADH_yield': nadh_yield,
            'ATP_from_NADH': atp_from_nadh,
            'Net_ATP_yield': net_atp,
            'Intermediate': 'Dihydroxyacetone phosphate (DHAP)'
        }

