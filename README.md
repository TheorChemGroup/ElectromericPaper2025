# Resurrecting The Electromeric Effect: Dynamic Polarizability Controls Selectivity

This repository contains all the raw computational data and a Python script to reproduce the plots presented in the associated manuscript. 
It supports reproducibility and transparency by executing the data analysis from the raw data which led us to our conclusions.

This project depends on GoodVibes2 software. Please ensure GoodVibes2 is installed and accessible in your system’s PATH before running the Python script.
For detailed installation instructions, please refer to the official GoodVibes2 repository: https://github.com/TheorChemGroup/GoodVibes2

## SI_archive Structure

- `reproduce.py`: Main Python script to process data and generate plots.
- `data/`: Subfolder containing prepared raw CSV data files and XYZ file with located stationary points.
-`Energy_profiles/`: Subfolder containing ORCA 5.0.4 `.log` files of 3', TS1, IM, TS2, 4' states calculations.
- `E_HU/`: Subfolder containing ORCA 5.0.4 `.log` files of simplified_IM structure calculations with orbital data printed for molecular orbital visualization.
- `Hammett_constants`: Subfolder containing ORCA 5.0.4 `.log` files necessary for Hammett constants calculation.

## Usage

1. Clone or download this repository:
    ```
    git clone https://github.com/TheorChemGroup/ElectromericPaper2025
    ```
2. Navigate to the `SI_archive` folder:
    ```
    cd SI_archive
    ```
3. Run the script to reproduce plots:
    ```
    python reproduce.py
    ```

## Output files

After running `reproduce.py`, the following files will be generated in the working directory:

- `Fig_*.csv`: Data files corresponding to figures presented in the manuscript. 
- `Fig_S21.csv`: Structural parameters (substituent, charge, multiplicity, number of imaginary frequencies) 
and energies of HOMO, LUMO, LUMO+1 orbitals. HOMO-LUMO Gap column is included for convenience.
- `Hammett_data_summary.csv`: File with structure parameters (substituent, charge, multiplicity, number of imaginary frequencies) 
and electronic energies. qhG energies are present in the last column and used for Hammett Constants calculation.
- `Structural_summary_data.csv`: Structural parameters, electronic and qhG energies.
Relative qhG energies are present in the last column and used for energy profiles analysis.
- `SI_structures_trj.xyz`: XYZ coordinates of all located stationary points

## Computational Details

For detailed computational methodology and additional information, please refer to the Supplementary Information.

## Contact

For questions or issues, please open an issue on GitHub or contact xdzhiblavi@gmail.com (Khadizha Dzhiblavi).




