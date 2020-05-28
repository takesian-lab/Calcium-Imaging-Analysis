# Calcium-Imaging-Analysis

Authors: Anne Takesian, Carolyn Sweeney, Maryse Thomas, and Wisam Reid

##### This is a work in progress


The stim_protocol code is:
- noiseburst      = 1
- ReceptiveField  = 2
- FM sweep        = 3
- widefield       = 4
- SAM             = 5
- SAM freq        = 6
- behavior        = 7
- behavior        = 8
- random h20      = 9
- noiseburst_ITI  = 10
- random air      = 11

-----

# Running the code
## Command line

1. Compile blocks
- Make Info excel spreadsheet outside of MATLAB. Look at the example spreadsheet and "Data Structure for Info" PDF within this folder for help.
- Update paths at the top of compile_blocks_from_info
- Run compile_blocks_from_info
- Visualize the output of a single block with the function visualize_block

2. Analyze
- You can use either the same Info spreadsheet as before or create a dedicated spreadsheet for each type of analysis you want to do.
- Update paths at the top of anlaysis_2P_v3
- Run analysis_2P_v3

# Code Structure

-----

# Development

Search **TODO** in the source code files


| TODO: | Code Location              | Call Location              | Task                                                             |
|:-----:|:--------------------------:|:--------------------------:|:----------------------------------------------------------------:|
| 1.    | [example.m]      | [main.m] | add a feature                                 |


## Miscellaneous Notes

# Dependencies
