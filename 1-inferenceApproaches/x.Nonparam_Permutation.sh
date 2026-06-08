#!/bin/bash

#######################################
# Choose options below (1=run)
#######################################
DO_Randomise=1
DO_func2anat=1
#######################################

parent_dir="" 
task="AUDIO"
version="ConsOrth_concat_sm" 
model="ICA" 
thresh_opt="vox_based" #vox_based, tfce

subject=""

for mask_location in wholeBrain thalamus_thr50 cerebellum brainstem; do
  if [ $mask_location == 'wholeBrain' ] ; then
    mask="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_run1/desc-optcom_bold_mask"
  elif [ $mask_location == 'thalamus_thr50' ] ; then
    mask_name="harvardoxford-subcortical_prob_CombinedThalamus_thr50_CORRECT_dilateSitekMGN_bin"
    mask="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin"
  elif [ $mask_location == 'brainstem' ] ; then
    mask_name="harvardoxford-subcortical_Brainstem_thr25_extend-IC_CORRECT_bin"
    mask="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin"
  elif [ $mask_location == 'cerebellum' ] ; then
    mask_name="mni_prob_Cerebellum_p50_2mm_bin"
    mask="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin"
  fi

output_folder="${parent_dir}/derivatives/${subject}/${session}/func/${task}"
contrastName="randomise_Audio-0_${mask_location}_${model}_${version}"
output_dir=${output_folder}/${contrastName}

if [ ! -d ${output_dir} ]
then
  mkdir ${output_dir}
fi

#########################
####### Randomise ########
#########################

if [ "${DO_Randomise}" -eq 1 ]
then
  echo "Running Randomise"
  echo "*****************"

# TASK REGRESSORS 
for cond1 in Audio; do 
cond2="0" #for paired 2-sample t-test cond1-cond2 (0 otherwise)

perm=500 # 5000 (for final testing) - maximum number of permutations / = 2^datasets
runMerge="fslmerge -t ${output_dir}/allbetamaps.nii.gz"

# ## Audio
if [ ! -f ${output_dir}/allbetamaps.nii.gz ] ; then
  for run in run1 run2 run3 run4
  do
    if [ ! -f ${mask}.nii.gz ] ; then
      # Transform mask to functional space
      ./x.PreProc_Transform_nonlin.sh ${parent_dir}/masks/${mask_name} \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.bet/${subject}_AudioRun1-SBREF_1_dc_bet_ero \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_stand2func_warp \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm \
        stand2func
        
      fslmaths ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func.nii.gz -bin ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin.nii.gz # binarize mask
    fi

    if [ "${cond2}" -eq 0 ]
    then
      ### 1-sample t-test
      # combine all beta parameter maps
      bcoef="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_${task}_${run}_TaskReg_bcoef.nii.gz" 
      runMerge="${runMerge} ${bcoef}"
    fi
  done

  eval ${runMerge}
else
  echo "File allbetamaps already exists! - Running randomise"
fi

if [ $thresh_opt == 'tfce' ]; then # run randomise with TFCE
  randomise -i ${output_dir}/allbetamaps.nii.gz -o ${output_dir}/${contrastName} \
    -m ${mask} -d ${output_dir}/design.mat -t ${output_dir}/design.con  -e ${output_dir}/design.grp \
    -T -v 5 -n ${perm} --glm_output
  # remove -v 5 (variance smoothing with sigma 5mm) if sample size > 20
  # -T does the following: threshold-free cluster enhancement; image remains voxel-wise and "un-clustered"

elif [ $thresh_opt == 'vox_based' ]; then   # run randomise with voxel-based thresholding
  randomise -i ${output_dir}/allbetamaps.nii.gz -o ${output_dir}/${contrastName} \
    -m ${mask} -d ${output_dir}/design.mat -t ${output_dir}/design.con  -e ${output_dir}/design.grp \
    -x -v 5 -n ${perm} --glm_output
  # remove -v 5 (variance smoothing with sigma 5mm) if sample size > 20
  # -x does the following: uses the null distribution of the max voxel-wise test statistic

fi

done
fi

#########################
####### Func2Anat ########
#########################

if [ "${DO_func2anat}" -eq 1 ]
then
  echo "Running func2anat"
  echo "*****************"

  input_file="${output_folder}/${contrastName}/randomise_Audio-0_${mask_location}_${model}_${version}_glm_cope"
  if [ ! -f ${input_file}_func2anat.nii.gz ] ; then
        # Transform input file to anatomical space
        ./x.PreProc_Transform_lin.sh ${input_file} \
          ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
          ${output_folder}/${contrastName} \
          func2anat
  fi

  if [ $thresh_opt == 'tfce' ]; then
    tstat_file="${output_folder}/${contrastName}/randomise_Audio-0_${mask_location}_${model}_${version}_tfce_corrp_tstat1"
  elif [ $thresh_opt == 'vox_based' ]; then
    tstat_file="${output_folder}/${contrastName}/randomise_Audio-0_${mask_location}_${model}_${version}_vox_corrp_tstat1"
  fi

  if [ ! -f ${tstat_file}_func2anat.nii.gz ] ; then
        # Transform input file to anatomical space
        ./x.PreProc_Transform_lin.sh ${tstat_file} \
          ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
          ${output_folder}/${contrastName} \
          func2anat
  fi

fi

done #end mask_location loop