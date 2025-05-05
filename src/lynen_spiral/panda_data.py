import xml.etree.ElementTree as ET
import pandas as pd

# Parse the XML file
tree = ET.parse('hmdb_metabolites.xml')
root = tree.getroot()

# Namespace dictionary (update with the correct namespace if needed)
ns = {'hmdb': 'http://www.hmdb.ca'}

# List to store fatty acid data
fatty_acids = []

# Iterate through each metabolite entry
for metabolite in root.findall('hmdb:metabolite', ns):
    # Extract the name
    name = metabolite.find('hmdb:name', ns).text if metabolite.find('hmdb:name', ns) is not None else None
    
    # Extract the SMILES string
    smiles = metabolite.find('hmdb:smiles', ns).text if metabolite.find('hmdb:smiles', ns) is not None else None
    
    # Extract the superclass
    superclass = metabolite.find('hmdb:super_class', ns).text if metabolite.find('hmdb:super_class', ns) is not None else None
    
    # Check if the superclass indicates a fatty acid
    if superclass and 'Fatty Acids' in superclass:
        fatty_acids.append({'Name': name, 'SMILES': smiles})

# Create a pandas DataFrame
df = pd.DataFrame(fatty_acids)

# Display the DataFrame
print(df.head())
