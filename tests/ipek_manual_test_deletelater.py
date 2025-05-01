from rdkit import Chem
from rdkit.Chem import AllChem

# Thioesterized palmitic acid
mol = Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)SC")

# Try dehydrogenation
rxn_d = AllChem.ReactionFromSmarts("[C:1]-[CH2:2]-C(=O)-S >> [C:1]=[C:2]-C(=O)-S") 
products_d = rxn_d.RunReactants((mol,))
for i, ps in enumerate(products_d):
    d_product = ps[0]
    print(f"Product {i}: {Chem.MolToSmiles(d_product)}")

# Sanitize and prepare the dehydrogenation product
Chem.SanitizeMol(d_product)

# Try hydration
rxn_h = AllChem.ReactionFromSmarts("[C:1](=O)[C:2]=[C:3] >> [C:1](=O)[C:2]([OH])[C:3]")
products_h = rxn_h.RunReactants((d_product,))
for i, ps in enumerate(products_h):
    h_product = ps[0]
    print(f"\nProduct {i}: {Chem.MolToSmiles(h_product)}")

# Sanitize and prepare the hydration product
Chem.SanitizeMol(h_product)

# More specific oxidation pattern
# Ensuring that the hydroxyl group is adjacent to the carbonyl group
rxn_o = AllChem.ReactionFromSmarts("[C:1]-[C:2]([OH:3])-[C:4](=O)S >> [C:1](=O)-[C:2]-[C:4](=O)S")
products_o = rxn_o.RunReactants((h_product,))
for i, ps in enumerate(products_o):
    o_product = ps[0]
    print(f"\nProduct {i}: {Chem.MolToSmiles(o_product)}")

Chem.SanitizeMol(o_product)

# Try thiolysis
rxn_t = AllChem.ReactionFromSmarts("[C:1](=[O:2])[C:3][C:4](=[O:5])S >> [C:1](=[O:2])O")
products_t = rxn_t.RunReactants((o_product,))
for i, ps in enumerate(products_t):
    t_product = ps[i]
    print(f"\nProduct {i}: {Chem.MolToSmiles(ps[0])}")

