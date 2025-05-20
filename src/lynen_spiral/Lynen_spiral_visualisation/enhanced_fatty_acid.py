import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from enum import Enum
from itertools import zip_longest

import sys
import os
sys.path.append(os.path.abspath(r"C:\Users\ipeki\git\Lynen-s_Spiral\Lynen-s_Spiral\src\lynen_spiral"))
from fatty_acid import FattyAcidMetabolism, FattyAcidType #type: ignore


# SMARTS patterns for each reaction step

SMARTS_REACTIONS = {
    # 1. Dehydrogenation: Remove two hydrogens to form a double bond between alpha and beta carbon
    "dehydrogenation": "[C:1]-[CH2:2]-C(=O)-S >> [C:1]=[C:2]-C(=O)-S",

    # OR

    # Δ²–Δ³ enoyl‑CoA isomerase‑assisted dehydrogenation (simplified as "alternative dehydrogenation")
    # Occurs when a pre‑existing Cβ=Cγ double bond has to be shifted.
    "alternative_dehydrogenation": "[C:1]=[C:2]-[C:3]-C(=O)-S >> [C:1]-[C:2]=[C:3]-C(=O)-S",

    # 2. Hydration: Add OH to beta-carbon, H to alpha-carbon
    "hydration": "[C:1]=[C:2]~[C:3](=O)-S >> [C:1]-[C:2]([OH])-[C:3](=O)-S",

    # 3. Oxidation: Oxidize beta-carbon alcohol (OH) to ketone (=O)
    "oxidation": "[C:1]-[C:2]([OH:3])-[C:4](=O)S >> [C:1](=O)-[C:2]-[C:4](=O)S",

    # 4. Thiolysis: Cleave between alpha and beta carbon; attach CoA-S group
    # Note: Representing CoA simply as 'S' attached to C
    "thiolysis": "[C:1](=[O:2])[C:3][C:4](=[O:5])S >> [C:1](=[O:2])S.[C:3][C:4](=[O:5])S"
    
}


# EnhancedFattyAcidMetabolism inherits all the attributes and methods of the FattyAcidMetabolism class, 
# while also having the ability to introduce its own unique attributes and methods or override existing ones.

class EnhancedFattyAcidMetabolism(FattyAcidMetabolism):
    """Enhanced version of FattyAcidMetabolism with reaction SMARTS implementation."""
    
    def __init__(self, input_value):
        """
        Initialize the object using a SMILES string

        Parameters:
            input_value (str): Input representing the molecule, which can be a SMILES string, 
                           a common fatty acid name, or a delta notation string.

        Notes:
            Uses the "_process_input" method to convert the input into a valid SMILES string.
        """

        smiles = self._process_input(input_value)
        super().__init__(smiles)

        self.reaction_steps = []
        self.reaction_results = []
        self.reaction_descriptions = []
        self.cycle_log = []
        self.atp_yield = 0

    def _process_input(self, input_value):

        """
        Convert the given input into a valid SMILES string.

        Attempts to parse the input_value as a SMILES string. If successful, returns the input as-is.
        If parsing fails or the input is invalid, raises a ValueError.

        Parameters:
            input_value (str): The input representing a molecule; expected to be a SMILES string.

        Returns:
            str: A valid SMILES string corresponding to the input.

        Raises:
            ValueError: If the input cannot be parsed into a valid SMILES string or if
                        the format is unrecognized.
        """

        try:
            # If input is a valid SMILES string
            mol = Chem.MolFromSmiles(input_value)
            if mol is not None:
                return input_value

        except Exception as e:
            raise ValueError(f"Error processing input '{input_value}': {str(e)}")

        # Input is unrecognized
        raise ValueError(
            f"Unrecognized input format: {input_value}. Use SMILES."
        )

    # ------------------ Methods to Hande Cis Double Bonds ------------------

    def get_double_bond_info(self, mol):
        
        """
        Identify all double bonds in the molecule and retrieve their properties.

        For each double bond, gathers:
        - The atom indices involved in the bond.
        - The position of the double bond within the molecule (using a helper method).
        - The stereochemistry of the bond ('cis' or 'trans').
        - The indices of neighboring atoms adjacent to the double bond atoms (used for E/Z determination).

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object to analyze.

        Returns:
            list of dict: A list where each dictionary describes a double bond with keys:
                - 'atoms' (tuple): Indices of the two atoms forming the double bond.
                - 'position' (int): Position of the double bond in the molecule.
                - 'stereo' (str): Either 'cis' or 'trans' stereochemistry.
                - 'neighbors' (list of list): Neighbor atom indices for each bond atom, 
                                              excluding the bonded partner.

        Notes:
            Uses the "_get_double_bond_position" method to determine bond position.
        """
    
        double_bonds = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Get the atoms involved in the double bond
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
            
                # Determine if it's cis or trans
                stereo = bond.GetStereo()
                is_cis = (stereo == Chem.BondStereo.STEREOCIS)
            
                # Get neighbor atoms for E/Z determination
                neighbors = []
                for atom in [begin_atom, end_atom]:
                    neighbors.append([n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() not in [begin_atom.GetIdx(), end_atom.GetIdx()]])
            
                position = self._get_double_bond_position(mol, begin_atom.GetIdx())
                bond_info = {
                    'atoms': (begin_atom.GetIdx(), end_atom.GetIdx()),
                    'position': position,
                    'stereo': 'cis' if is_cis else 'trans',
                    'neighbors': neighbors
                }
                double_bonds.append(bond_info)

        return double_bonds

    def _get_double_bond_position(self, mol, atom_idx):

        """
        Determine the position of a double bond relative to the thioester end (carbonyl carbon) 
        of the fatty acid.

        Starting from the given atom index (one atom of the double bond), this method traverses 
        along the carbon chain towards the thioester end, counting carbons until it reaches the 
        terminal point.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object.
            atom_idx (int): The index of one atom in the double bond whose position is being 
                            determined.

        Returns:
            int: The position of the double bond counted as the number of carbons from the 
                 thioester carbon (starting at 1).

        Notes:
            Assumes a linear carbon chain starting from the thioester carbon.
        """

        # Traverse from the thioester carbon (C=O-S), starting from position 0
        position = 0
        atom = mol.GetAtomWithIdx(atom_idx)
    
        # Count atoms from the carbonyl carbon
        while atom.GetAtomicNum() == 6:  # While it's a carbon
            position += 1
            # Get the neighbor that's part of the main chain
            neighbors = [n for n in atom.GetNeighbors() 
                        if n.GetAtomicNum() == 6 and n.GetIdx() != atom_idx]
            if not neighbors:
                break
            atom_idx = atom.GetIdx()
            atom = neighbors[0]
    
        return position

    def handle_cis_double_bonds(self, mol):

        """
        Detect and convert any Δ³-cis double bonds in the molecule to the trans configuration.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object to process.

        Returns:
            rdkit.Chem.Mol: The updated molecule with Δ³-cis double bonds converted to trans.
                            If no such bonds are found, returns the original molecule unchanged.

        Notes:
            - Only the Δ³-cis double bonds are modified; other cis or trans bonds remain unchanged.
            - Uses the "get_double_bond_info" method to identify double bonds and their stereochemistry.
            - Uses the "_convert_cis_to_trans" method to perform the conversion.
        """

        double_bonds = self.get_double_bond_info(mol)
    
        for bond_info in double_bonds:
            # Get the actual RDKit bond
            bond = mol.GetBondBetweenAtoms(*bond_info['atoms'])
            
            # Only act on Δ³-cis regardless of even/odd
            if bond_info['stereo'] == 'cis' and bond_info['position'] == 3:
                return self._convert_cis_to_trans(mol, bond_info['atoms'])
    
        return mol

    def _convert_cis_to_trans(self, mol, atom_indices):

        """
        Convert a specified cis double bond in the molecule to trans configuration.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object to modify.
            atom_indices (tuple of int): Pair of atom indices defining the double bond.

        Returns:
            rdkit.Chem.Mol: The updated molecule with the specified double bond converted to trans.
        """

        # Find the bond between cis double-bond atoms
        bond = mol.GetBondBetweenAtoms(atom_indices[0], atom_indices[1])
        if bond:
            # Set the bond stereochemistry to trans
            bond.SetStereo(Chem.BondStereo.STEREOTRANS)
            # Update the molecule
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol)
        
            # Log this modification
            self.reaction_steps.append("Cis-Trans Isomerization")
            self.reaction_descriptions.append(
                f"Converted cis double bond (atoms {atom_indices}) to trans configuration"
            )
        return mol
    
    # ------------------ Methods to Handle Double Bonds Between Cα-Cβ or Cβ-Cγ ------------------
    
    def find_carbonyl_carbon(self, mol):

        """
        Locate the carbonyl carbon atom (C=O) specifically in the thioester functional group (C=O–S).

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object to analyze.

        Returns:
            int or None: The atom index of the carbonyl carbon if found; otherwise, None.
        """

        # Iterate through all atoms in the molecule
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "S":
                # Check for a neighboring carbon with a double bond to oxygen (C=O)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "C":
                        # Check if this carbon has a double bond to oxygen
                        for carbonyl_neighbor in neighbor.GetNeighbors():
                            if carbonyl_neighbor.GetSymbol() == "O" and mol.GetBondBetweenAtoms(neighbor.GetIdx(), carbonyl_neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                return neighbor.GetIdx()  # Return the index of the carbonyl carbon
        return None  # Return None if no carbonyl carbon is found
    
    def get_alpha_beta_gamma_carbons(self, mol):

        """
        Identify and return the indices of the alpha (Cα), beta (Cβ), and gamma (Cγ) carbon atoms
        relative to a carbonyl carbon in the given molecule.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object.

        Returns:
            tuple: A tuple (c_alpha, c_beta, c_gamma) of atom indices.
               If any of the positions can't be found, the respective entry will be None.
               For example, (index_of_Cα, index_of_Cβ, None) if Cγ is missing.
        
        Notes:
            Uses the "find_carbonyl_carbon" method.
        """

        carbonyl_c = self.find_carbonyl_carbon(mol)
        if carbonyl_c is None:
            return None, None, None

        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_c)

        # Find Cα: carbon directly bonded to the carbonyl carbon
        c_alpha = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Carbon
                c_alpha = nbr.GetIdx()
                break

        if c_alpha is None:
            return None, None, None

        # Find Cβ: carbon bonded to Cα, but not the carbonyl carbon
        c_beta = None
        for nbr in mol.GetAtomWithIdx(c_alpha).GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_c:
                c_beta = nbr.GetIdx()
                break

        if c_beta is None:
            return c_alpha, None, None

        # Find Cγ: carbon bonded to Cβ, not Cα
        c_gamma = None
        for nbr in mol.GetAtomWithIdx(c_beta).GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c_alpha:
                c_gamma = nbr.GetIdx()
                break

        return c_alpha, c_beta, c_gamma

    def has_alpha_beta_double_bond(self, mol):

        """
        Check if there is a double bond between the alpha (Cα) and beta (Cβ) carbon atoms
        relative to the carbonyl carbon in the molecule.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object.

        Returns:
            bool: True if an α-β double bond exists, False otherwise.

        Notes: 
            Uses the "get_alpha_beta_gamma_carbons" method.
        """

        c_alpha, c_beta, c_gamma = self.get_alpha_beta_gamma_carbons(mol)
        if c_alpha is None or c_beta is None:
            return False

        # Check α-β double bond
        bond_ab = mol.GetBondBetweenAtoms(c_alpha, c_beta)
        if bond_ab and bond_ab.GetBondType() == Chem.BondType.DOUBLE:
            return True

        return False
    
    def has_beta_gamma_double_bond(self, mol):

        """
        Check if there is a double bond between the beta (Cβ) and gamma (Cγ) carbon atoms
        relative to the carbonyl carbon in the molecule.

        Parameters:
            mol (rdkit.Chem.Mol): The RDKit molecule object.

        Returns:
            bool: True if a β-γ double bond exists, False otherwise.
        
            Notes: 
            Uses the "get_alpha_beta_gamma_carbons" method
        """
        
        c_alpha, c_beta, c_gamma = self.get_alpha_beta_gamma_carbons(mol)
        if c_alpha is None or c_beta is None:
            return False

        if c_gamma is not None:
            bond_bg = mol.GetBondBetweenAtoms(c_beta, c_gamma)
            if bond_bg and bond_bg.GetBondType() == Chem.BondType.DOUBLE:
                return True
        
        return False

    # ------------------ Methods to Simulate Activation and Transport ------------------

    def activate_fatty_acid(self):

        """
        Activate a fatty acid by converting it into an acyl-CoA derivative.
        This simulates the biological activation step via acyl-CoA synthetase,
        consuming the equivalent of 2 ATP molecules.

        The method applies a reaction SMARTS to simulate the attachment of CoA 
        (represented here as a sulfur atom, as RDKit doesn't recognize ) to the carboxylic acid group.

        Returns:
            rdkit.Chem.Mol: The RDKit molecule object representing the activated fatty acid (acyl-CoA).

        Raises:
            RuntimeError: If the reaction fails or produces no product.
        """

        try:
            # SMARTS pattern for the activation reaction, representing
            # the conversion of a carboxylic acid to an acyl-CoA
            activation_rxn = AllChem.ReactionFromSmarts("[C:1](=[O:2])[O:3] >> [C:1](=[O:2])[S:3]")

            if not activation_rxn:
                raise ValueError("Invalid activation reaction SMARTS pattern")
        
            # Prepare reactant
            reactant = self.molecule
        
            # Run reaction
            products = activation_rxn.RunReactants((reactant,))
        
            if not products or len(products) == 0:
                raise RuntimeError("Activation reaction failed to produce any products")
        
            # Get the product molecule
            activated_mol = products[0][0]

            # Sanitize the product molecule
            activated_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(activated_mol)

            # Store reaction information
            self.reaction_steps.append("Activation")
            self.reaction_results.append(activated_mol)
            self.reaction_descriptions.append("Conversion to acyl-CoA via acyl-CoA synthetase (2 ATP cost)")
        
            return activated_mol
        
        except Exception as e:
            # Provide detailed error information
            raise RuntimeError(f"Failed ot activate fatty acid. Error: {str(e)}\n"
                               f"Input Molecule: {Chem.MolToSmiles(self.molecule, canonical=True)}")
    
    
    def carnitine_shuttle(self, acyl_coa_mol):

        """
        Simulate the mitochondrial transport of an acyl-CoA molecule via the carnitine shuttle.

        In cellular metabolism, long-chain fatty acids must be transported into the mitochondrial 
        matrix for β-oxidation. This transport is facilitated by the carnitine shuttle, which 
        involves a series of enzymatic steps:

            1. Acyl-CoA reacts with carnitine to form acylcarnitine and free CoA (via CPT1 on the outer mitochondrial membrane).
            2. Acylcarnitine is translocated across the inner mitochondrial membrane.
            3. Acylcarnitine reacts with CoA to regenerate acyl-CoA and release carnitine (via CPT2 on the inner membrane).

        For simplicity and clarity, these intermediate steps are **not explicitly modeled** in this method.
        Instead, the method assumes successful transport and proceeds directly to β-oxidation.

        Parameters:
            acyl_coa_mol (rdkit.Chem.Mol): The RDKit molecule representing the fatty acyl-CoA to be transported.

        Returns:
            rdkit.Chem.Mol: The same molecule object, now considered to be inside the mitochondrial matrix.
        """
        
        # Store reaction information
        self.reaction_steps.append("Transport")
        self.reaction_results.append(acyl_coa_mol)  # No molecular change, just location
        self.reaction_descriptions.append("Transport into mitochondria via carnitine shuttle")
        
        return acyl_coa_mol
    
    # ------------------------------- Methods to Run Reactions -------------------------------

    def run_reaction(self, reaction_smarts, molecule):

        """
        Apply a SMARTS-based chemical reaction to an RDKit molecule.

        This method compiles a reaction from a SMARTS pattern and applies it to the given molecule.
        It is designed to support both single-product and multi-product reactions, such as those
        encountered in fatty acid metabolism (e.g., cleavage of β-ketoacyl-CoA into acetyl-CoA 
        and a shortened acyl-CoA).

        The method handles reaction sanitization, basic error checking, and product retrieval.
        If two products are produced (as in β-oxidation), it logs both structures and returns
        the first one for further processing.

        This method applied in the beta_oxidation_cycle method.

        Parameters:
            reaction_smarts (str): A SMARTS reaction string describing the transformation.
            molecule (rdkit.Chem.Mol): The RDKit molecule to which the reaction will be applied.

        Returns:
            rdkit.Chem.Mol or None: The primary product of the reaction if successful, or None 
            if the reaction fails or produces no products.

        Raises:
            RuntimeError: If the reaction SMARTS is invalid or the transformation fails during execution.
        """

        try:
            # Compile the reaction from SMARTS
            reaction = AllChem.ReactionFromSmarts(reaction_smarts)
            if not reaction:
                raise ValueError(f"Invalid SMARTS reaction pattern: {reaction_smarts}")
            
            # Run the reaction
            products = reaction.RunReactants((molecule,))

            if not products or len(products) == 0:
                print(f"Warning: Reaction produced no products. Reaction SMARTS: {reaction_smarts}")
                return None  # Or handle appropriately

            if len(products) > 0 and len(products[0]) > 0:
                # Check if the reaction produced more than one product
                if len(products) > 1:
                    product_1 = products[0][0]  # First product (acetyl-CoA)
                    product_2 = products[0][1]  # Second product (shortened-CoA)

                    try:
                        Chem.SanitizeMol(product_1)
                        molecule.UpdatePropertyCache(strict=False) # Sanitize before the reaction
                
                    except Exception as e:
                        print(f"Error sanitizing molecule before reaction: {e}")
                        return None  # Or handle appropriately
            
                    return product_1
                
                else:
                    molecule = products[0][0] # Access the first molecule of the first product set
                    
                try:
                    Chem.SanitizeMol(molecule)
                    molecule.UpdatePropertyCache(strict=False) # Sanitize before the reaction
                
                except Exception as e:
                    print(f"Error sanitizing molecule before reaction: {e}")
                    return None  # Or handle appropriately
                    
                return molecule
            
            else:
                print(f"No products generated. Reaction SMARTS: {reaction_smarts}")
                return None 

        except Exception as e:
            # provide detailked error information
            raise RuntimeError(f"Failed to run reaction. Error: {str(e)}\n"
                               f"Reaction SMARTS: {reaction_smarts}\n"
                               f"Molecule: {Chem.MolToSmiles(molecule, canonical=True)}")
    

    # Example of one beta-oxidation cycle with reaction SMARTS
    def beta_oxidation_cycle(self, acyl_coa_mol):

        """
        Perform one complete β‑oxidation cycle on an acyl‑CoA molecule.

        The routine automatically decides between the standard
        Cα–Cβ dehydrogenation and an alternative dehydrogenation
        (Δ²–Δ³ isomerase‑assisted) whenever a pre‑existing Cβ=Cγ bond is
        detected.  After dehydrogenation, the canonical hydration,
        oxidation, and thiolytic cleavage steps are executed.

        This method is applied in the run_complete_oxidation method.

        Reaction sequence implemented
        -----------------------------
        1. Dehydrogenation 
            • Enzyme : acyl‑CoA dehydrogenase (FAD → FADH₂, +1.5 ATP eq)  
            • If a Cβ=Cγ double bond already exists, the alternative_dehydrogenation
            SMARTS is used instead. (No ATP yield in this case.)

        2. Hydration
            • Enzyme : enoyl‑CoA hydratase (crotonase)  
            • No direct ATP yield.

        3. Oxidation  
            • Enzyme : 3‑hydroxyacyl‑CoA dehydrogenase (NAD⁺ → NADH, +2.5 ATP eq).

        4. Thiolysis
            • Enzyme : β‑ketoacyl‑CoA thiolase  
            • Releases acetyl‑CoA and shortens the chain by two carbons.

        Parameters
            acyl_coa_mol : rdkit.Chem.Mol
                The RDKit molecule representing the input fatty acyl‑CoA inside
                the mitochondrial matrix.

        Returns
            rdkit.Chem.Mol
                The shortened acyl‑CoA produced at the end of the cycle
                (ready for another round of β‑oxidation).

        Raises
            RuntimeError
                Propagated from `run_reaction` if a SMARTS pattern is invalid
                or a reaction fails.

        Notes:
            Uses the following methods:
            - "handle_cis_double_bonds"
            - "has_beta_gamma_double_bond"
            - "has_alpha_beta_double_bond"
            - "run_reaction"
        """
    
        # Check for and handle any problematic double bonds
        acyl_coa_mol = self.handle_cis_double_bonds(acyl_coa_mol)

        # Define reaction steps
        beta_reaction_steps = []

        # Check for existing double bonds and decide on the reaction steps
        if self.has_beta_gamma_double_bond(acyl_coa_mol):

            beta_reaction_steps.append(("Isomerase‑assisted dehydrogenation", 
                                        SMARTS_REACTIONS["alternative_dehydrogenation"], 1.5))
            self.reaction_steps.append("Isomerase‑assisted dehydrogenation")
            self.reaction_results.append(acyl_coa_mol)
            self.reaction_descriptions.append("Isomerase‑assisted dehydrogenation due to existing double bond at Δ² or Δ³")
        
        elif not self.has_alpha_beta_double_bond(acyl_coa_mol):

            beta_reaction_steps.append(("Dehydrogenation", SMARTS_REACTIONS["dehydrogenation"], 1.5)) # ATP yield for FADH2
            
        # Define a list of the canonical reaction steps using the SMARTS patterns
        beta_reaction_steps.extend([
            ("Hydration", SMARTS_REACTIONS["hydration"], 0),        # No ATP yield
            ("Oxidation", SMARTS_REACTIONS["oxidation"], 2.5),      # ATP yield for NADH
            ("Thiolysis", SMARTS_REACTIONS["thiolysis"], 0),        # No ATP yield
        ])
    
        # Store intermediate products and descriptions
        products = [acyl_coa_mol]  # Initialize with input molecule
        cycle_steps = []  # Start new cycle log

        for step_name, reaction_smarts, atp_yield in beta_reaction_steps:

            # Run each reaction
            current_mol = products[-1]  # Get the most recent product

            try: 
                Chem.SanitizeMol(current_mol)  # Ensure the molecule is valid
                current_mol.UpdatePropertyCache(strict=False)  # Ensure implicit valence is calculated

            except Exception as e:
                print(f"Error sanitizing molecule at step {step_name}: {e}")
                for atom in current_mol.GetAtoms():
                    print(f"Atom {atom.GetIdx()} - {atom.GetSymbol()} valence: {atom.GetTotalValence()}")
                continue

            # Run the reaction
            result = self.run_reaction(reaction_smarts, current_mol)

            try:
                Chem.SanitizeMol(result)
                result.UpdatePropertyCache()
                
            except Exception as e:
                print(f"Error sanitizing result at step {step_name}: {e}")
                continue

            if result is None:
                print(f"Invalid result at step {step_name}")
                continue

            # Store reaction information
            self.reaction_steps.append(step_name)
            self.reaction_results.append(result)
            self.reaction_descriptions.append(f"{step_name} (SMARTS: {reaction_smarts})")

            # Log ATP yield
            self.atp_yield += atp_yield

            # Append new product and cycle step information
            products.append(result)

            cycle_steps.append({
                "step": step_name,
                "input": Chem.MolToSmiles(current_mol),
                "output": Chem.MolToSmiles(result),
                "smarts": reaction_smarts,
                "atp_yield": atp_yield
                })

        # Log the cycle information
        self.cycle_log.append({
            "cycle_number": len(self.cycle_log) + 1,
            "steps": cycle_steps
            })

        return products[-1]  # Return the final product of the cycle
    
    def run_complete_oxidation(self):
        
        """
        Core method of the EnhancedFattyAcid Class.
        Run the complete beta-oxidation pathway for a fatty acid molecule.

        This method simulates the full process of fatty acid degradation:
    
        1. Activation of the fatty acid to acyl-CoA (2 ATP equivalent cost).
        2. Transport of the acyl-CoA into the mitochondria via the carnitine shuttle.
        3. Iterative beta-oxidation cycles, each removing a 2-carbon acetyl-CoA unit and 
            generating ATP from FADH₂ and NADH.
        4. Handles:
            - Even-chain fatty acids: Fully converted to acetyl-CoA.
            - Odd-chain fatty acids: Ends with propionyl-CoA.
            - Edge case: Molecule too short to oxidize (failsafe condition).

        The method tracks:
            - All intermediate steps and transformations.
            - Total ATP yield from FADH₂ and NADH (β-oxidation ATP yield).
            - Reaction logs including SMILES of intermediate molecules.

        Returns:
            dict: A summary of the oxidation process with keys:
                - 'final_products': {
                        'acetyl_coa_count': (int) number of acetyl-CoA units formed,
                        'propionyl_coa': (bool) True if a 3-carbon propionyl-CoA was produced
                        }
                - 'reaction_steps': (list) names of steps performed
                - 'reaction_results': (list) RDKit Mol objects for each step
                - 'reaction_descriptions': (list) human-readable reaction explanations
                - 'total_atp_yield': (float) net ATP yield (FADH₂ + NADH equivalents only)

        Raises:
            RuntimeError: If activation, transport, or a beta-oxidation step fails.

        Notes:
            - This does not simulate TCA cycle oxidation of acetyl-CoA or gluconeogenesis
              from propionyl-CoA.
            - Propionyl-CoA is only noted as present.
            - ATP yield includes only direct cofactors from beta-oxidation.
            - Uses the following methods:
                - "activate_fatty_acid"
                - "carnitine_shuttle"
                - "_get_carbon_chain_length"
                - "beta_oxidation_cycle"
                - "calculate_atp_yield"
        """

        prev_length = None

        # Clear previous results
        self.reaction_steps = []
        self.reaction_results = []
        self.reaction_descriptions = []
        self.reaction_cycles = [] 
        self.cycle_log = []
        self.atp_yield = 0  # Initialize ATP yield

        # Initial step: Activation
        current_mol = self.activate_fatty_acid()

        # Transport step: Carnitine shuttle
        current_mol = self.carnitine_shuttle(current_mol)

        # Collect products
        acetyl_coa_count = 0
        propionyl_coa = None
        cycle = 0

        while True:
            chain_length = self._get_carbon_chain_length(current_mol)

            # Handle a chain length that isn't decreasing
            # (Should not happen under normal circumstances)
            if chain_length == prev_length:
                print("Chain length not decreasing — breaking to avoid infinite loop.")
                break

            prev_length = chain_length

            # Handle final steps
            if chain_length == 2:
                # Final step for even chain: acetyl-CoA
                self.reaction_cycles.append(f"Final step: 2-carbon fragment to acetyl-CoA")
                self.reaction_descriptions.append("Last 2-carbon unit converted to acetyl-CoA")
                acetyl_coa_count += 1
                break
            elif chain_length == 3:
                # Final step for odd-chain: propionyl-CoA
                self.reaction_cycles.append(f"Final step: 3-carbon fragment to propionyl-CoA")
                self.reaction_descriptions.append("Last 3-carbon unit converted to propionyl-CoA")
                propionyl_coa = current_mol
                break
            elif chain_length < 2: 
                # Handle edge case: chain too short to process 
                # (Should not happen under normal circumstances)
                self.reaction_cycles.append("Error: Chain too short to process")
                break

            # Run a regular beta-oxidation cycle
            cycle += 1
            self.reaction_cycles.append(f"Cycle {cycle}")
            self.reaction_descriptions.append(f"Starting beta-oxidation cycle {cycle}")
            current_mol = self.beta_oxidation_cycle(current_mol)
            acetyl_coa_count += 1

        # Final energy calculation
        self.atp_yield += self.calculate_atp_yield()['total_ATP']

        return {
            'final_products': {
                'acetyl_coa_count': acetyl_coa_count,
                'propionyl_coa': propionyl_coa is not None
            },
            'reaction_steps': self.reaction_steps,
            'reaction_results': self.reaction_results,
            'reaction_descriptions': self.reaction_descriptions,
            'total_atp_yield': self.atp_yield
        }

    # ----------------------------- Method to Prepare Data for Visualization -----------------------------

    def prepare_data_for_visualization(self):

        """
        Prepare a detailed, modular data structure suitable for visualization.

        Processes the reaction steps, descriptions, and resulting molecules to
        generate a list of step dictionaries containing metadata for each step,
        including molecular formula and canonical SMILES representation.

        Returns:
            dict: A dictionary with the following structure:
                {
                    "metadata": {
                        "total_steps": int,       # Number of reaction steps
                        "total_atp_yield": float, # Total ATP yield from the process
                        "total_cycles": int       # Number of cycles completed
                        },

                    "steps": [                   # List of dicts for each reaction step
                        {
                            "index": int,        # Zero-based step index
                            "step_number": int,  # One-based step number
                            "name": str,         # Name of the reaction step
                            "description": str,  # Description of the step
                            "smiles": str,       # Canonical SMILES string of the molecule
                            "formula": str       # Molecular formula
                        },
                        ...
                    ],
                    "cycles": list              # Raw cycle log data from the oxidation process
                }

        Notes:
            - If a molecule cannot be sanitized, the formula and SMILES will be set to "Invalid".
            - The "cycles" entry provides detailed cycle info useful for downstream visualization.

        """
        
        steps = []

        for idx, (step, description, molecule) in enumerate(
            zip(self.reaction_steps, self.reaction_descriptions, self.reaction_results)
        ):
            try:
                Chem.SanitizeMol(molecule)
                formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)
                smiles = Chem.MolToSmiles(molecule, canonical=True)
            except Exception:
                formula = "Invalid"
                smiles = "Invalid"

            steps.append({ # Each element in the list is a dictionary
                "index": int(idx),
                "step_number": idx + 1,
                "name": step,
                "description": description,
                "smiles": smiles,
                "formula": formula,
            })

            # Example of the structure of the returned data:

            # steps = [] → A list of dictionaries.
            # steps[0] → The first step’s dictionary.
            # steps[0]["name"] → The name of that first step, like "dehydrogenation"

        return {
            "metadata": {
                "total_steps": len(steps),
                "total_atp_yield": self.atp_yield,
                "total_cycles": len(self.cycle_log),
            },
            "steps": steps,
            "cycles": self.cycle_log
        }