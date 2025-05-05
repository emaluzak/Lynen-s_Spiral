from rdkit import Chem

import sys
import os

sys.path.append(os.path.abspath(r"C:\Users\ipeki\git\Lynen-s_Spiral\Lynen-s_Spiral\src\lynen_spiral"))

from enhanced_fatty_acid import EnhancedFattyAcidMetabolism #type: ignore

def manual_test_visualization_data():
    # Create molecule for palmitic acid
    #palmitic_acid = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)O")

    # Initialize the metabolism object
    pathway = EnhancedFattyAcidMetabolism("CCCCCCCCCCCCCCCC(=O)O")

    # Run your simulation (you may have your own wrapper method)
    #acyl_coa = pathway.activate_fatty_acid()
    pathway.run_complete_oxidation()

    # Get the visualization data
    viz_data = pathway.prepare_data_for_visualization()

    # Print results
    from pprint import pprint
    pprint(viz_data)

if __name__ == "__main__":
    manual_test_visualization_data()
