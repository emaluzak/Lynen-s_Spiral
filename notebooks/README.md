# Lynen's Spiral: Fatty Acid Beta-Oxidation Visualizer
**Practical Programming in Chemistry Project**  


**Lynen’s Spiral: Fatty Acid Beta-Oxidation Viusalizer** is an interactive Streamlit-based web application developed as part of the Practical Programming in Chemisty course. It is designed to facilitate a deeper understanding of the biochemical process of fatty acid beta-oxidation. By integrating molecular visualizations and stepwise reaction breakdowns, this tool aims to make complex metabolic pathways more accessible to students, educators, and researchers.

---

## Features

### 1. Homepage – Fatty Acid Selection

- Provides introductory context on fatty acids and the beta-oxidation pathway.
- Allows users to select from a library of predefined fatty acids or to construct custom molecules by specifying chain length, degree of unsaturation, and cis/trans configurations.
- Offers 2D and 3D molecular visualizations of the selected fatty acid.
- Selection unlocks access to subsequent simulation and visualization modules (Lynen Spiral and Mechanisms).

### 2. Lynen Spiral – 3D Visualization of Beta-Oxidation Cycles

- Presents a dynamic 3D visualization of the cyclic steps involved in fatty acid beta-oxidation.
- Users can navigate through each oxidation cycle via a cycle slider.
- A step slider enables in-depth exploration of individual biochemical steps within a given cycle.
- Each step includes a corresponding 2D molecular structure of the metabolic intermediate.

### 3. Mechanisms – Comprehensive Cycle Analysis

- Displays a detailed breakdown of all biochemical reactions for the selected fatty acid.
- Provides structural representations of reactants and products at each stage.
- Includes ATP yield calculations for complete fatty acid oxidation.

---

## Technical Architecture

- **Frontend and Visualization:** Streamlit  
- **3D Molecular Rendering:** Plotly  
- **2D Molecular Representations:** RDKit  
- **Backend Logic:** Python-based simulation of enzymatic processes and reactions using SMILES and reaction SMARTS

---

## Educational Objectives

This application is intended to support:

- A step-by-step conceptualization of fatty acid beta-oxidation.
- Interactive visualization of complex metabolic pathways in both 2D and 3D formats.
- Comparative analysis of how molecular characteristics (e.g., chain length, unsaturation) influence oxidation mechanisms and energetic output.

---

## Getting Started

### Prerequisites

- Python 3.8 or higher  
- `pip` (Python package manager)

### Installation

```bash
# Clone the repository
git clone https://github.com/emaluzak/Lynen-s_Spiral.git
cd Lynen-s_Spiral

# Create and activate a virtual environment (recommended)
python -m venv venv

# On Windows
venv\Scripts\activate

# On macOS/Linux
source venv/bin/activate

# Install required dependencies
pip install -r requirements.txt
```

---

## Dependencies

The following Python libraries are required:

- streamlit  
- plotly  
- rdkit  
- numpy  
- pandas  
- matplotlib 

---

## Running the Application

To launch the application, run the following command:

```bash
streamlit run Home.py
```

The application will open in your default web browser at:  
[http://localhost:8501](http://localhost:8501)

---

## Usage Recommendations

- Begin by selecting a fatty acid from the homepage or creating a custom molecule.
- Utilize the navigation sidebar to explore different modules.
- For optimal visualization, use a modern browser with WebGL support.
- For the best viewing experience, enable the light theme and "Wide Mode" via the Streamlit settings menu (accessible via the three-dot icon in the top-right corner).

---

## License

MIT License

---

## Acknowledgements

This project was developed as part of the Practical Programming in Chemistry curriculum to enhance biochemical education through interactive tools.  
Special thanks to faculty and peers who provided guidance during development.

## Authors

- [Ema Lužáková](https://github.com/emaluzak)
- [Ariane Fienga Schultheis](https://github.com/arianeschultheis)
- [Ipek Begüm Incesu](https://github.com/ipekinc)
