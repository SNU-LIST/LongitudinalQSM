# LongitudinalQSM
Joint reconstruction framework for longitudinal QSM that improves inter-scan consistency while preserving sensitivity to true susceptibility changes.

# Longitudinal QSM (MATLAB)

MATLAB implementation of a longitudinal QSM reconstruction framework based on joint multi-time-point estimation, designed to improve inter-scan consistency while preserving sensitivity to true susceptibility changes.

## Environment setup
This code is implemented in MATLAB.

The reconstruction is developed based on the MEDI QSM framework.  
Please install the MEDI toolbox prior to running the code:

MEDI toolbox: https://github.com/liu-lab/MEDI_toolbox

Make sure all required paths are properly added in MATLAB.

## Usage
Run the main script to perform longitudinal QSM reconstruction.

## Input
The input data should include:
- Local field maps for each time point
- Corresponding brain masks
- Magnitude images (for structural regularization)
- B0 direction
- voxel size

Optional:
- Refined brain mask excluding CSF and major vessels
- 
All inputs are expected to be preprocessed (e.g., phase unwrapping, background field removal, and registration across time points).

## License
© Jiye Kim  
@ LIST, Seoul National University  

This code is provided for academic and research purposes only.  
Please cite or acknowledge the original author when using this work.
