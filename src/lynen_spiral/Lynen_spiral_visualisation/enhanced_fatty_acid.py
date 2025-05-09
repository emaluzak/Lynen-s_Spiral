import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from enum import Enum
from itertools import zip_longest

# Import Ema's code (saved in a file called fatty_acid.py)
import sys
import os
sys.path.append(os.path.abspath(r"C:\Users\ipeki\git\Lynen-s_Spiral\Lynen-s_Spiral\src\lynen_spiral"))
from fatty_acid import FattyAcidMetabolism, FattyAcidType #type: ignore

# Common fatty acid names to SMILES dictionary
FATTY_ACIDS = {
    "palmitic acid": "CCCCCCCCCCCCCCCC(=O)O",    # C16:0
    "stearic acid": "CCCCCCCCCCCCCCCCCC(=O)O",   # C18:0
    "oleic acid": "CCCCCCCCC=CCCCCCCC(=O)O",     # C18:1(9)
    "linoleic acid": "CCCCCC=CCC=CCCCCCC(=O)O",  # C18:2(9,12)
    "alpha-linolenic acid": "CCC=CCC=CCC=CCCCCCC(=O)O",  # C18:3(9,12,15)
    "arachidonic acid": "CCCC=CCC=CCC=CCC=CCCCC(=O)O"    # C20:4(5,8,11,14)
}

# SMARTS patterns for each reaction step

SMARTS_REACTIONS = {
    # 1. Dehydrogenation: Remove two hydrogens to form a double bond between alpha and beta carbon
    "dehydrogenation": "[C:1]-[CH2:2]-C(=O)-S >> [C:1]=[C:2]-C(=O)-S",

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
        Initialize with SMILES, fatty acid name, or notation like C18:2(9,12).
        """
        smiles = self._process_input(input_value)
        print(f"Processed SMILES: {smiles}")  # Debugging output
        super().__init__(smiles)

        self.reaction_steps = []
        self.reaction_results = []
        self.reaction_descriptions = []
        self.cycle_log = []
        self.atp_yield = 0

    def _process_input(self, input_value):
        """Process input into a valid SMILES string."""
        try:
            # If input is a valid SMILES string
            mol = Chem.MolFromSmiles(input_value)
            if mol is not None:
                return input_value
        
            # If input matches a common fatty acid name (case-insensitive)
            input_value_lower = input_value.lower()
            if input_value_lower in FATTY_ACIDS:
                return FATTY_ACIDS[input_value_lower]

            # If input is in CX:Y(pos1,pos2,...) notation
            if input_value.startswith("C") and ":" in input_value:
                return self._notation_to_smiles(input_value)

        except Exception as e:
            raise ValueError(f"Error processing input '{input_value}': {str(e)}")

        # Input is unrecognized
        raise ValueError(
            f"Unrecognized input format: {input_value}. Use SMILES, common name, or CX:Y notation."
        )


    def parse_fatty_acid_notation(self, notation):
        """
        Parse and validate fatty acid notation like 'C18:2(9,12)'.
        Returns the number of carbons and a list of double bond positions.

        :param notation: String, fatty acid notation
        :return: (num_carbons, bond_positions)
        """
        try:
            # Split into parts
            parts = notation.split(':')
            num_carbons = int(parts[0][1:])  # Extract the number of carbons (e.g., "C18" -> 18)

            # Handle saturated fatty acids (no double bonds)
            if len(parts) == 1 or '(' not in parts[1]:
                return num_carbons, 0, []

            # Unsaturated: Extract double bond info
            double_bond_info = parts[1]
            double_bond_count = int(double_bond_info.split('(')[0])  # Number of double bonds
            bond_positions = [int(pos) for pos in double_bond_info.split('(')[1][:-1].split(',')]

            # Validate positions
            if double_bond_count != len(bond_positions):
                raise ValueError("Mismatch between bond count and positions.")
            if any(pos < 1 or pos >= num_carbons for pos in bond_positions):
                raise ValueError("Invalid bond position detected.")

            return num_carbons, double_bond_count, bond_positions

        except Exception as e:
            raise ValueError(f"Invalid fatty acid notation '{notation}': {str(e)}")
    

    def construct_backbone(self, num_carbons, double_bond_count, bond_positions):
        """
        Construct the carbon backbone of the fatty acid with double bonds.

        :param num_carbons: Number of carbons in the chain
        :param bond_positions: List of positions for double bonds
        :return: String representing the backbone SMILES (no carboxyl group yet)
        """

        bond_positions_2 = [x - double_bond_count for x in bond_positions]
        backbone = []
        for i in range(1, num_carbons + 1 - double_bond_count):
            if i in bond_positions_2:
                backbone.append("C=")
            else: 
                backbone.append("C")  # Append a single-bonded carbon

        return ''.join(backbone)


    def add_carboxyl_group(self, backbone):
        """
        Add the carboxyl group to the backbone.

        :param backbone: String, the carbon backbone
        :return: Complete SMILES string with carboxyl group
        """
        reversed_backbone = backbone[::-1]
        return reversed_backbone + "(=O)O"

    def _notation_to_smiles(self, notation):
        """
        Convert CX:Y(pos1,pos2,...) notation to SMILES.

        :param notation: String, fatty acid notation
        :return: SMILES string for the fatty acid molecule
        """
        try:
            # Parse input
            num_carbons, double_bond_count, bond_positions = self.parse_fatty_acid_notation(notation)

            # Construct backbone
            backbone = self.construct_backbone(num_carbons, double_bond_count, bond_positions)

            # Add carboxyl group
            return self.add_carboxyl_group(backbone)

        except Exception as e:
            raise ValueError(f"Error converting notation '{notation}': {str(e)}")

    # ------------------ Methods to Hande Cis Double Bonds ------------------

    def get_double_bond_info(self, mol):
        """Identify all double bonds and their properties (position, cis/trans)."""
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

                print(f"Atom {bond_info['atoms'][0]}-{bond_info['atoms'][1]} | "
                      f"Position: {bond_info['position']} | "
                      f"Stereo: {bond_info['stereo']}")

        return double_bonds

    def _get_double_bond_position(self, mol, atom_idx):
        """Determine the position of a double bond relative to the thioester end."""
        # Traverse from the thioester carbon (C=O-S) to find position
        # This assumes the thioester is at position 1
        position = 0
        atom = mol.GetAtomWithIdx(atom_idx)
    
        # Simple approach - count atoms from the carbonyl carbon
        # You may need to adjust this based on your exact molecular representation
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
        """Convert Δ³-cis bonds to trans."""
        double_bonds = self.get_double_bond_info(mol)
    
        for bond_info in double_bonds:
            # Get the actual RDKit bond
            bond = mol.GetBondBetweenAtoms(*bond_info['atoms'])
        
            # Debug print to verify everything
            print(f"\nBond between atoms {bond_info['atoms']}:")
            print(f"From dict: position={bond_info['position']}, stereo={bond_info['stereo']}")
            print(f"From RDKit: stereo={bond.GetStereo()} (cis={bond.GetStereo() == Chem.BondStereo.STEREOCIS})")
            
            # Only act on Δ³-cis regardless of even/odd
            if bond_info['stereo'] == 'cis' and bond_info['position'] == 3:
                print(f"Converting Δ³-cis to Δ²-trans at atoms {bond_info['atoms']}")
                return self._convert_cis_to_trans(mol, bond_info['atoms'])
    
        return mol

    def _convert_cis_to_trans(self, mol, atom_indices):
        """Convert a specific cis double bond to trans configuration."""
        # Find the bond between these atoms
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
        Find the carbonyl carbon (C=O) in the thioester group (C=O–S).
        Returns the index of the carbonyl carbon atom if found, otherwise None.
        """
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "S":
                # Check for a neighboring carbon with a double bond to oxygen (C=O)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "C":
                        # Check if this carbon has a double bond to oxygen
                        for carbonyl_neighbor in neighbor.GetNeighbors():
                            if carbonyl_neighbor.GetSymbol() == "O" and mol.GetBondBetweenAtoms(neighbor.GetIdx(), carbonyl_neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                print(f"Found carbonyl carbon at index {neighbor.GetIdx()} (C=O-S) with sulfur at index {atom.GetIdx()}")
                                return neighbor.GetIdx()  # Return the index of the carbonyl carbon
        return None  # Return None if no carbonyl carbon is found
    
    def get_alpha_beta_gamma_carbons(self, mol):
        carbonyl_c = self.find_carbonyl_carbon(mol)
        if carbonyl_c is None:
            print("No carbonyl carbon found.")
            return None, None, None

        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_c)

        # Find Cα: carbon directly bonded to the carbonyl carbon
        c_alpha = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Carbon
                c_alpha = nbr.GetIdx()
                print(f"Found Cα at index {c_alpha}")
                break

        if c_alpha is None:
            print("No Cα found.")
            return None, None, None

        # Find Cβ: carbon bonded to Cα, but not the carbonyl carbon
        c_beta = None
        for nbr in mol.GetAtomWithIdx(c_alpha).GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_c:
                c_beta = nbr.GetIdx()
                print(f"Found Cβ at index {c_beta}")
                break

        if c_beta is None:
            print("No Cβ found.")
            return c_alpha, None, None

        # Find Cγ: carbon bonded to Cβ, not Cα
        c_gamma = None
        for nbr in mol.GetAtomWithIdx(c_beta).GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c_alpha:
                c_gamma = nbr.GetIdx()
                print(f"Found Cγ at index {c_gamma}")
                break

        return c_alpha, c_beta, c_gamma

    def has_alpha_beta_double_bond(self, mol):
        c_alpha, c_beta, c_gamma = self.get_alpha_beta_gamma_carbons(mol)
        if c_alpha is None or c_beta is None:
            return False

        # Check α-β double bond
        bond_ab = mol.GetBondBetweenAtoms(c_alpha, c_beta)
        if bond_ab and bond_ab.GetBondType() == Chem.BondType.DOUBLE:
            print(f"Found double bond between Cα (index {c_alpha}) and Cβ (index {c_beta}).")
            return True

        print(f"No double bond found between Cα-Cβ.")
        return False
    
    def has_beta_gamma_double_bond(self, mol):
        c_alpha, c_beta, c_gamma = self.get_alpha_beta_gamma_carbons(mol)
        if c_alpha is None or c_beta is None:
            return False

        if c_gamma is not None:
            bond_bg = mol.GetBondBetweenAtoms(c_beta, c_gamma)
            if bond_bg and bond_bg.GetBondType() == Chem.BondType.DOUBLE:
                print(f"Found double bond between Cβ (index {c_beta}) and Cγ (index {c_gamma}).")
                return True
        
        print(f"No double bond found between Cβ-Cγ.")
        return False

    # Reaction SMARTS implementation for activation (acyl-CoA synthetase) TO MODIFY !!!
    def activate_fatty_acid(self):
        """
        Activate fatty acid by converting to acyl-CoA.
        This process consumes 2 ATP equivalent.
        """

        try:
            # Simplified SMARTS pattern to fix activation reaction
           
            # Replace this line:
            # activation_rxn = AllChem.ReactionFromSmarts("[C:1](=[O:2])[O:3][H] >> [C:1](=[O:2])[S:4][CoA]")

            # With a simpler version:
            activation_rxn = AllChem.ReactionFromSmarts("[C:1](=[O:2])[O:3] >> [C:1](=[O:2])[S:3]C")


            if not activation_rxn:
                raise ValueError("Invalid activation reaction SMARTS pattern")
        
            # Prepare reactant
            reactant = self.molecule

            # Debugging: Print reactant molecule to check its structure
            print("Reactant Molecule:", Chem.MolToSmiles(reactant, canonical=True))
        
            # Run reaction
            products = activation_rxn.RunReactants((reactant,))
        
            if not products or len(products) == 0:
                raise RuntimeError("Activation reaction failed to produce any products")
        
            # Get the product molecule
            activated_mol = products[0][0]

            # Debugging: Print product molecule to check if it's formed
            print("Activated Molecule:", Chem.MolToSmiles(activated_mol, canonical=True))

            # Fix before adding hydrogens
            activated_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(activated_mol)
        
            # Add hydrogens for better visualization
            #activated_mol = Chem.AddHs(activated_mol)
        
            # Store reaction information
            self.reaction_steps.append("Activation")
            self.reaction_results.append(activated_mol)
            self.reaction_descriptions.append("Conversion to acyl-CoA via acyl-CoA synthetase (2 ATP cost)")
            print(f"Activation step: {self.reaction_steps}")
        
            return activated_mol
        
        except Exception as e:
            # Provide detailed error information
            raise RuntimeError(f"Failed ot activate fatty acid. Error: {str(e)}\n"
                               f"Input Molecule: {Chem.MolToSmiles(self.molecule, canonical=True)}")
    
    # Implement carnitine shuttle
    def carnitine_shuttle(self, acyl_coa_mol):
        """
        Simulate transport of acyl-CoA into mitochondria via carnitine shuttle.
        This is a multi-step process:
        1. Acyl-CoA + Carnitine → Acylcarnitine + CoA (via CPT1)
        2. Transport across membrane
        3. Acylcarnitine + CoA → Acyl-CoA + Carnitine (via CPT2)
        """
        # For simplicity, we'll simulate this without actual reaction SMARTS
        # In a real implementation, you'd model each step
        
        # Store reaction information
        self.reaction_steps.append("Transport")
        self.reaction_results.append(acyl_coa_mol)  # No molecular change, just location
        self.reaction_descriptions.append("Transport into mitochondria via carnitine shuttle")
        print(f"After transportation step: {self.reaction_steps}")
        
        return acyl_coa_mol
    
    def run_reaction(self, reaction_smarts, molecule):
        """Run a SMARTS-based reaction."""
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
                    print(f"Product 1 (Acetyl-CoA): {Chem.MolToSmiles(product_1)}")
                    print(f"Product 2 (Shortened-CoA): {Chem.MolToSmiles(product_2)}")

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
                return None  # Handle appropriately, e.g., return None or raise an error

        except Exception as e:
            # provide detailked error information
            raise RuntimeError(f"Failed to run reaction. Error: {str(e)}\n"
                               f"Reaction SMARTS: {reaction_smarts}\n"
                               f"Molecule: {Chem.MolToSmiles(molecule, canonical=True)}")
    

    # Example of one beta-oxidation cycle with reaction SMARTS
    def beta_oxidation_cycle(self, acyl_coa_mol):
        """
        Perform one complete beta-oxidation cycle with cis/trans handling:

        1. Dehydrogenation, 
            Formation of a double bond between C_alpha and C_beta. 
            (FAD → FADH₂)
            Product: Beta-dehydroacyl-CoA
            Catalyzed by: Flavin dehydrogenase (acyl-CoA dehydrogenase)

        2. Hydration
            Addition of H2O to the double bond
            Product: beta-hydroxyacyl-CoA
            Catalyzed by: Lyase (enoyl-CoA hydratase, also called crotonase)

        3. Oxidation 
            (aka second dehydrogenation)
            (NAD⁺ → NADH)
            Product: beta-oxoacyl-CoA (beta-ketoacyl-CoA)
            Catalyzed by: Pyridine dehydrogenase (3-hydroxyacyl-CoA dehydrogenase)

        4. Thiolysis
            The beta-oxoacyl-CoA is highly unstable as a thioester and undergoes thiolytic cleavage.
            Catalyzed by: Acyltransferase (beta-oxothiolase)
            End result: Acetyl-CoA is released, and the fatty acid chain is shortened by two carbon atoms.
        """
    
        # First check for and handle any problematic double bonds
        acyl_coa_mol = self.handle_cis_double_bonds(acyl_coa_mol)
        print("\n--- After Isomerization ---")
        print("Molecule:", Chem.MolToSmiles(acyl_coa_mol))

        # Define reaction steps
        beta_reaction_steps = []


        if self.has_beta_gamma_double_bond(acyl_coa_mol):

            beta_reaction_steps.append(("Alternative Dehydrogenation", SMARTS_REACTIONS["alternative_dehydrogenation"], 2))
            print("Alternative Dehydrogenation: Existing double bond between Cβ and Cγ")
            self.reaction_steps.append("Alternative Dehydrogenation")
            self.reaction_results.append(acyl_coa_mol)
            self.reaction_descriptions.append("Alternative dehydrogenation due to existing double bond at Δ² or Δ³")
        
        elif not self.has_alpha_beta_double_bond(acyl_coa_mol):

            beta_reaction_steps.append(("Dehydrogenation", SMARTS_REACTIONS["dehydrogenation"], 2))

        else:
            print("Double bond between Cα and Cβ, skipping Hydrogenation.")
            
        # Define a list of the reaction steps using the SMARTS patterns
        beta_reaction_steps.extend([
            ("Hydration", SMARTS_REACTIONS["hydration"], 0),                # No ATP yield
            ("Oxidation", SMARTS_REACTIONS["oxidation"], 3),                # ATP yield for NADH
            ("Thiolysis", SMARTS_REACTIONS["thiolysis"], 0),                # No ATP yield
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

                print(f"Before {step_name}: {Chem.MolToSmiles(current_mol)}")

            except Exception as e:
                print(f"Error sanitizing molecule at step {step_name}: {e}")
                for atom in current_mol.GetAtoms():
                    print(f"Atom {atom.GetIdx()} - {atom.GetSymbol()} valence: {atom.GetTotalValence()}")
                continue

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

            print(f"After {step_name}: {Chem.MolToSmiles(result)}")

            # If multiple products (e.g., thiolysis), choose the longer fragment
            #if isinstance(result, tuple):
                # For thiolysis, pick the longer fatty acid fragment (avoid acetyl-CoA)
                #result = max(result, key=lambda mol: mol.GetNumAtoms())

            # Store reaction information
            self.reaction_steps.append(step_name)
            self.reaction_results.append(result)
            self.reaction_descriptions.append(f"{step_name} (SMARTS: {reaction_smarts})")
            print(f"\nAfter {step_name} step: {self.reaction_steps}")

            # Log ATP yield
            self.atp_yield += atp_yield

            # Append new product
            products.append(result)

            cycle_steps.append({
                "step": step_name,
                "input": Chem.MolToSmiles(current_mol),
                "output": Chem.MolToSmiles(result),
                "smarts": reaction_smarts,
                "atp_yield": atp_yield
            })

        self.cycle_log.append({
            "cycle_number": len(self.cycle_log) + 1,
            "steps": cycle_steps
     })

        return products[-1]  # Return the final product of the cycle
    
    def run_complete_oxidation(self):
        """
        Run the complete beta-oxidation process from start to finish.
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
        print(Chem.MolToSmiles(current_mol))  # Log the molecule state

        # Transport step: Carnitine shuttle
        current_mol = self.carnitine_shuttle(current_mol)
        print(Chem.MolToSmiles(current_mol))  # Log the molecule state

        # Collect products
        acetyl_coa_count = 0
        propionyl_coa = None
        cycle = 0

        while True:
            chain_length = self._get_carbon_chain_length(current_mol)
            print(f"Cycle {cycle}: Chain length = {chain_length}")

            if chain_length == prev_length:
                print("Chain length not decreasing — breaking to avoid infinite loop.")
                break

            prev_length = chain_length

            # Handle final steps
            if chain_length == 2:
                # Final acetyl-CoA
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

    
    def prepare_data_for_visualization(self):
        """Prepare a detailed, modular data structure for visualization."""
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

            steps.append({ # Appending dictionaries for each step!
                "index": int(idx),
                "step_number": idx + 1,
                "name": step,
                "description": description,
                "smiles": smiles,
                "formula": formula,
            })

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
            "cycles": self.cycle_log,  # <- added raw cycle log directly
        }
    

    def visualize_reaction_sequence(self):
        """Visualize the β-oxidation reaction steps."""
        steps = []  # List to hold molecule images for each step

        # Automatically generate the reaction_cycles list
        num_non_cycle_steps = 2  # "activation" and "transport" steps
        total_steps = num_non_cycle_steps + len(self.reaction_steps)

        # Prepare the visualization for each step
        for idx, (step, description, molecule) in enumerate(
            zip(
                self.reaction_steps,
                self.reaction_descriptions,
                self.reaction_results,
            )
        ):
            # Generate a molecule image for this step
            try:
                mol_image = Draw.MolToImage(molecule, size=(300, 300))
            except Exception as e:
                print(f"Error generating image for step {idx+1}: {e}")
                continue

            # Store step data for visualization
            steps.append({
                "step_number": idx + 1, 
                "step_name": step,
                "description": description,
                "molecule_image": mol_image
            })

        # Visualize all the steps (this could be displayed or saved depending on your visualization platform)
        fig, axs = plt.subplots(1, len(steps), figsize=(5 * len(steps), 5))
        for ax, step in zip(axs, steps):
            ax.imshow(step["molecule_image"])
            ax.set_title(f"Step {step['step_number']}: {step['step_name']}")
            ax.axis('off')

        plt.tight_layout()
        plt.show()
