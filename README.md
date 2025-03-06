# Protein Quality Check Program  

## Overview  
The **Protein Quality Check Program** provides functionality to **automatically download** a protein from the **Protein Data Bank (PDB)** and analyze its **polypeptide regions** for potential **low-quality** areas.  

The analysis is based on several key criteria, including:  
- **X-ray data resolution**  
- **Refinement metrics** (R-factor & R-free)  
- **B-factors (temperature factors)**  
- **Atomic occupancy**  
- **Model geometry** (Ramachandran plot)  

## Implementation Details  

The program is implemented using an **object-oriented approach** with the `Quality` class.  
### `Quality` Class  
#### Constructor:  
```python
def __init__(self, structure, path):

Initializes the class with structure and path as inputs.
Key Methods:
poly_detail(self)

Connects to the polypeptide builder and retrieves polypeptide chains and their lengths.
get_resolution(self)

Extracts the resolution of the structure from structure.header["resolution"].
get_r_factor(self)

Retrieves the R-factor by converting mmCIF to a dictionary using:
python
Copy
Edit
build_dict = MMCIF2Dict(path)[8]
get_data(self)

Displays overall protein data, including:
Polychain number
Residue number & name
Atoms
B-factor
Occupancy
Torsion angles
Chain ID
cal_bad_b_factor_region(self) & cal_b_factor(self)

Identifies regions with high B-factors (greater than 80).
occupancy(self)

Checks for atoms with occupancy < 1.
rachamandran_outliers(self) & rachamandran_outliers_region(self)

Identifies amino residues in Ramachandran's disallowed region.
low_quality_region(self)

Combines results from B-factor, occupancy, and Ramachandran analysis to identify low-quality regions and provides details such as:
Chain ID
Polypeptide number
Residue ID
Residue name
Potential issue (e.g., high B-factor, low occupancy, or Ramachandran outlier)
give_score(self)

Calculates a quality score and classifies the structure quality as:
High Quality: Score = 5
Medium Quality: Score = 3
Low Quality: Score = 1

### Usage
git clone https://github.com/Tinsusu/structural-bioinformatics.git
