# Description

My solution to the first Budweiser problem: Data Conversion and the Preparation of a Purity Report for the given data. 

This is incomplete, and will be updated as we obtain both the descriptions and data to complete it.

It may also be updated as errors are found or instructions are updated.

# Usage

1. You can run either `python3 generate_key.py` or run all cells in the `generate_key.ipynb` to get the same output - both produce the `key.csv` file given the single input `StandardsData.xlsx` file. If you have a different input filename than `StandardsData.xlsx`, instead run `python3 generate_key.py InputStandardsFile.xlsx`, e.g. for `StandardsProfiles.xlsx`, run `python3 generate_key.py StandardsProfiles.xlsx`
2. Modify the key.csv file if you wish to add further LINE_NM -> PROFILE mappings.
3. You can then run either `python3 generate_report.py` or run all cells in the `generate_report.ipynb` to get the same final output - both produce the `report.csv` and `supplemental.csv` files given the input `GenoResults.csv` file and the `key.csv` file created by Step #1. If you have a different input filename than `GenoResults.csv`, instead run `python3 generate_key.py InputGenoFile.csv`, e.g. for `GenoMatrix.csv`, run `python3 generate_key.py GenoMatrix.csv`

Good luck, have fun,

-Blake Edwards / Dark Element
