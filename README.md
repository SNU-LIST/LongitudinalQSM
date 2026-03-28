
# Longitudinal QSM (MATLAB)

MATLAB implementation of a longitudinal QSM reconstruction framework based on joint multi-time-point estimation, designed to improve inter-scan consistency while preserving sensitivity to true susceptibility changes.

## Environment setup
This code is implemented in MATLAB.

The reconstruction is developed based on the MEDI QSM framework.  
Please install the MEDI toolbox prior to running the code:

MEDI toolbox: https://pre.weill.cornell.edu/mri/pages/qsm.html

Make sure all required paths are properly added in MATLAB.

## Usage
Run the main script to perform longitudinal QSM reconstruction of two time points.

## Input

### Required
- `RDF1`, `RDF2`  
  Local field maps in radians  
  (convert from Hz using `2π * delta_TE`)

- `Mask1`, `Mask2`  
  Brain masks

- `voxel_size`  
  Voxel size (mm)

- `delta_TE`  
  Echo spacing (s)

- `CF`  
  Center frequency (Hz)

- `B0_dir1`, `B0_dir2`  
  B0 direction vectors (unit vectors)
- `iMag1`, `iMag2`  
  Magnitude images for MEDI regularization

- `VesselMask`  
  Mask to exclude vessels from temporal sparsity

- `Mask_CSF`  
  CSF mask for susceptibility referencing

- `x0_init_ppm`, `dx_init_ppm`  
  Initial maps (ppm) (derived from MEDI)

All inputs are expected to be preprocessed (e.g., phase unwrapping, background field removal, and registration across time points).

## License
© Jiye Kim  
@ LIST, Seoul National University  

This code is provided for academic and research purposes only.  
Please cite or acknowledge the original author when using this work.
