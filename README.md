# PrecisionMapping_Statistics
Analysis code to accompany "Evaluating Approaches for Inference Testing of Whole-Brain Densely Sampled Single-Subject Task fMRI Data" paper

* Subject-level analysis code can be found on the following Github repository (https://github.com/BrightLab-ANVIL/WholeBrain_AuditoryMapping) on folder 1-subjectGLM.

## 1. Analysis for Inference Approaches

**Whole-Brain Analyses**
- x.WholeBrainVoxFDR.sh: Uses x.Subject_SigActivation.sh script to find voxels of significant positive activation with False Discovery Rate (FDR) correction

- x.WholeBrain_ClustFWE.sh: Example script to perform cluster-level inference testing with different thresholds

- x.Nonparam_Permutation.sh: Example script to perform non-parametric permutation inference testing with voxel_level (vox_based) or threshold free cluster enhancement (tfce) thresholding. Use mask_location = wholeBrain.

**Region-Specific Analyses**
- x.RegionSpecific_Top%tstat.sh: Uses DO_TopPer option to find voxels above a specified threshold (e.g. top 5% voxels) based on t-statistic value

- x.Nonparam_Permutation.sh: Example script to perform non-parametric permutation inference testing with voxel_level (vox_based) or threshold free cluster enhancement (tfce) thresholding. Use specific ROI for the mask_location loop (e.g. brainstem).
  
## 2. Stability Calculations

- Folds_COM.R: Example script to calculate the euclidean distance between the center of mass (COM) of an inference method's 40-minute data activation cluster and the COM of its respective 30-minute folds.

- Folds_DiceCoeff.R: Example script to calculate the dice coefficient score between an inference method's 40-minute data activation cluster and its respective 30-minute folds. 
  
