#!/bin/bash

#######################################
# Choose analysis options below (1=run)
#######################################
DO_SigActivation=1

parent_dir=""
mask_dir=""
whole_brain="MNI152_T1_2mm_brain"
task="AUDIO"
run_num=4 #Number of concatenated runs
ses_num=1 #If runs are from the same session ses_num=1; otherwise ses_num=2
fold_id=d #a, b, c, d
model="ICA"

## CHOOSE version
## With all 40 min of data => ConsOrth_concat${run_num}runs_sm_${ses_num}ses 
## With 30-minute folds => ConsOrth_concat${run_num}${fold_id}runs_sm_${ses_num}ses
version="ConsOrth_concat${run_num}runs_sm_${ses_num}ses" 
 

################################
################################
################################

for subject in sub-01 sub-02 sub-03
do

for session in ses-02
do

output_folder="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}"
space="func" # standard, func

echo "*********************"
echo "Processing ${subject}"
echo "*********************"

## CHOOSE threshold
for alpha in 05 01 001 # corresponding to 0.05, 0.01, and 0.001, respectively
do

if [ "${DO_SigActivation}" -eq 1 ]

then
  echo "****************************"
  echo "Running subject-level significant activation"
  echo "****************************"

  input_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat"
  brain_mask=${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_run1/desc-optcom_bold_mask.nii.gz
  matrix_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_matrix.1D"
  ROI="wholeBrain"

  sig_file="${input_file}_fdrp${alpha}_${ROI}"
  if [ ! -f ${sig_file}.nii.gz ] ; then
    ./x.Subject_SigActivation.sh ${input_file} ${brain_mask} ${matrix_file} ${alpha} ${ROI}
  fi
  
  if [ ! -f ${sig_file}_func2anat.nii.gz ] ; then
    # Transform output file (sig_file) to anatomical space
    ./x.PreProc_Transform_lin.sh ${sig_file} \
      ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
      ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
      ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
      func2anat
          
    fslmaths ${sig_file}_func2anat.nii.gz -bin ${sig_file}_func2anat_bin.nii.gz # binarize mask
  fi

fi

done # end alpha loop
done # end session loop
done # end subject loop