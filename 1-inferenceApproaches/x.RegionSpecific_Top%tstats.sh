#!/bin/bash

#######################################
# Choose analysis options below (1=run)
#######################################
DO_TopPer=1
DO_TopPer_transform=1

parent_dir=""
mask_dir=""
whole_brain="MNI152_T1_2mm_brain"
task="AUDIO"
run_num=4 #Number of concatenated runs
ses_num=1 
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
space="func"

echo "*********************"
echo "Processing ${subject}"
echo "*********************"

## CHOOSE threshold
for above in 95 90 # 90(top 10%), 95(top 5%)
do

if [ "${DO_TopPer}" -eq 1 ]

then
  echo "****************************"
  echo "Running Extract top %"
  echo "****************************"


  for ROI in thalamus_thr50 cerebellum brainstem; do
    if [ $ROI == 'thalamus_thr50' ] ; then
      mask_name="harvardoxford-subcortical_prob_CombinedThalamus_thr50_CORRECT_dilateSitekMGN_bin"
    elif [ $ROI == 'brainstem' ] ; then
      mask_name="harvardoxford-subcortical_Brainstem_thr25_extend-IC_CORRECT_bin"
    elif [ $ROI == 'cerebellum' ] ; then
      mask_name="mni_prob_Cerebellum_p50_2mm_bin"
    fi

  if [ $space == 'func' ] ; then
    mask="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin"

    if [ ! -f ${mask}.nii.gz ] ; then
      # Transform mask to functional space
      ./x.PreProc_Transform_nonlin.sh ${mask_dir}/${mask_name} \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.bet/${subject}_AudioRun1-SBREF_1_dc_bet_ero \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_stand2func_warp \
        ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm \
        stand2func

      fslmaths ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func.nii.gz -bin ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${mask_name}_stand2func_bin.nii.gz # binarize mask
    fi

    input_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat"
    output_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat_${ROI}_thrP${above}"

    if [ ! -f ${input_file}_${ROI}.nii.gz ] ; then
      3dcalc -a ${input_file}.nii.gz -b ${mask}.nii.gz -expr 'a*b' -prefix ${input_file}_${ROI}.nii.gz # extract mask area of tstat
    fi

    bcoef_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_bcoef"
    if [ ! -f ${bcoef_file}_${ROI}.nii.gz ] ; then
      3dcalc -a ${bcoef_file}.nii.gz -b ${mask}.nii.gz -expr 'a*b' -prefix ${bcoef_file}_${ROI}.nii.gz # extract mask area of bcoef
    fi
  
    fslmaths ${input_file}_${ROI}.nii.gz -thrP ${above} ${output_file}.nii.gz # find top X% of tstats within mask  
  fi
  done # end ROI loop
fi

if [ "${DO_TopPer_transform}" -eq 1 ]

then
  echo "****************************"
  echo "Running transform top %"
  echo "****************************"


  for ROI in thalamus_thr50 cerebellum brainstem; do 
    if [ $ROI == 'thalamus_thr50' ] ; then
      mask_name="harvardoxford-subcortical_prob_CombinedThalamus_thr50_CORRECT_dilateSitekMGN_bin"
    elif [ $ROI == 'brainstem' ] ; then
      mask_name="harvardoxford-subcortical_Brainstem_thr25_extend-IC_CORRECT_bin"
    elif [ $ROI == 'cerebellum' ] ; then
      mask_name="mni_prob_Cerebellum_p50_2mm_bin"
    fi
    

  above=95 #95(top 5%)

  bcoef_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_bcoef"
  if [ ! -f ${bcoef_file}_func2anat.nii.gz ] ; then
        # Transform output file (topPer) to anatomical space
        ./x.PreProc_Transform_lin.sh ${bcoef_file} \
          ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
          func2anat
  fi

  input_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat"
  if [ ! -f ${input_file}_func2anat.nii.gz ] ; then
        # Transform output file (topPer) to anatomical space
        ./x.PreProc_Transform_lin.sh ${input_file} \
          ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
          func2anat
  fi

  output_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat_${ROI}_thrP${above}"
  if [ ! -f ${output_file}_func2anat.nii.gz ] ; then
        # Transform output file (topPer) to anatomical space
        ./x.PreProc_Transform_lin.sh ${output_file} \
          ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
          ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
          func2anat
  fi

  if [ ! -f ${output_file}_func2stand.nii.gz ] ; then
        # Transform output file (topPer) to standard space
        ./x.PreProc_Transform_nonlin.sh ${output_file} \
            ${mask_dir}/${whole_brain} \
            ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2stand_warp \
            ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
            func2stand
  fi

  done # end ROI loop
fi

done ## end "above" threshold loop

done # end session loop
done # end subject loop