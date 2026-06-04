#!/bin/bash
#This function finds significant voxels on the subject level using outputs from the subject-level GLM
#Output files are created in the same directory as the input t-statistic map

#Check if the inputs are correct
if [ $# -ne 5 ]
then
  echo "Insufficient inputs"
  echo "Input 1 should be the t-statistic map (don't include '.nii.gz' in file name)"
  echo "Input 2 should be the brain mask (or ROI mask)"
  echo "Input 3 should be the matrix.1D file output from the subject-level GLM"
  echo "Input 4 should be the alpha level of significance (input 05 for alpha = 0.05)"
  echo "Input 5 should be the ROI name"
  exit
fi

input_tstat="${1}"
brain_mask="${2}"
matrix_file="${3}"
alpha="${4}"
ROI="${5}"

if [ ! -f "${input_tstat}_fdr${alpha}_${ROI}.nii.gz.txt" ]
then

  # find number of TRs and number of regressors from matrix.1D file
  matrix_rowcol=`1d_tool.py -show_rows_cols -infile ${matrix_file} -verb 0`
  arr=($matrix_rowcol)
  n_TRs=${arr[0]} # rows = number of TRs
  n_regressors=${arr[1]} # cols = number of regressors

  # calculate DOF
  ndof=$((${n_TRs}-${n_regressors}-1)) #N-k-1
  echo "ndof is $ndof"

  # find t-stat threshold for a < alpha
  tstat=$( cdf -p2t fitt 0.${alpha} ${ndof} )
  tstat=${tstat##* }
  echo "tstat is $tstat"

  # convert tstat to z score by FDR correction and threshold significant voxels
  3dFDR -input ${input_tstat}.nii.gz -mask ${brain_mask} -prefix ${input_tstat}_fdr_${ROI}.nii.gz # FDR correction
  fslmaths ${input_tstat}_fdr_${ROI}.nii.gz -thr ${tstat} ${input_tstat}_fdr${alpha}_${ROI}.nii.gz # threshold by tstat (includes pos and neg significant tstats)

  # find only significant POSITIVE voxels; remove this section if not needed
  fslmaths ${input_tstat}.nii.gz -thr 0 -bin ${input_tstat}_thrp_${ROI}.nii.gz # find where tstat is positive
    # to find where tstat is negative, use -uthr instead of -thr
  3dcalc -a ${input_tstat}_fdr${alpha}_${ROI}.nii.gz -b ${input_tstat}_thrp_${ROI}.nii.gz \
    -expr 'a*b' -prefix ${input_tstat}_fdrp${alpha}_${ROI}.nii.gz # find intersection of positive tstat and significant tstat maps --> only significant positive voxels

  # remove intermediate files
  rm ${input_tstat}_fdr_${ROI}.nii.gz ${input_tstat}_fdr${alpha}_${ROI}.nii.gz ${input_tstat}_thrp_${ROI}.nii.gz

else
  echo "** ALREADY RUN: tstat=${input_tstat} **"
fi
