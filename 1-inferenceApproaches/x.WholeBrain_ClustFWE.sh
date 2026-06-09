#!/bin/bash

parent_dir="" 
cond2="0"
run_num=3 #Choose based on how many concatenated runs are being used (i.e. 4 or 3 runs)
ses_num=1
fold_id=d #a, b, c, d
version_runs="ConsOrth_concat_sm" # Set name for individual runs
model="ICA" 
fxn="Cluster"
task="AUDIO"
athr="05" # alpha threshold for clustering

## CHOOSE version
## With all 40 min of data => ConsOrth_concat${run_num}runs_sm_${ses_num}ses 
## With 30-minute folds => ConsOrth_concat${run_num}${fold_id}runs_sm_${ses_num}ses
version="ConsOrth_concat${run_num}runs_sm_${ses_num}ses" 

for pthr in 005 001 # p-value threshold for clustering (005, 001)
do
  if [ ${pthr} == "005" ]; 
  then
    target_row=4 # 4 for p = 0.005, 6 for p = 0.001
    DO_acf=1
    DO_ClustSim=1
    DO_Clusterize=1
    DO_func2anat=1
  elif [ ${pthr} == "001" ]; 
  then
    target_row=6
    DO_acf=0
    DO_ClustSim=0
    DO_Clusterize=1
    DO_func2anat=1
  fi

for subject in sub-01
do
  echo "*********************"
  echo "Processing ${subject} with threshold p<${pthr}"
  echo "*********************"
  if [[ ${subject} == "sub-01" ]]; 
  then
    session="ses-04"
  fi

  brain_mask=${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_run1/desc-optcom_bold_mask.nii.gz # whole brain mask

for cond1 in Audio
do

  output_folder="${parent_dir}/derivatives/${subject}/${session}/func/${task}"
  contrastName="${fxn}_${cond1}-${cond2}_${model}_${version}"
  output_dir=${output_folder}/${contrastName}

  if [ ! -d ${output_dir} ]
  then
    mkdir ${output_dir}
  fi

  #########################
  ##### Calculate acf #####
  #########################

  if [ "${DO_acf}" -eq 1 ]
  then
    echo "Running acf"
    echo "*****************"

    # Calculate acf for each individual scan run using 3dFWHMx
    if [ ${run_num} == "4" ]; then
      echo "******************************************"
      echo "Running acf with 4 concatenated runs"
      echo "******************************************"
      for run in run1 run2 run3 run4 
      do
        run3dFWHMx="3dFWHMx -mask ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_${run}/desc-optcom_bold_mask.nii.gz"
        run3dFWHMx="${run3dFWHMx} -input ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version_runs}/${subject}_${task}_${run}_errts.nii.gz -acf ${output_dir}/${subject}_${task}_${run}_acf_emp.1D"
        run3dFWHMx="${run3dFWHMx} > ${output_dir}/${subject}_${task}_${run}_acf.1D" 
        eval ${run3dFWHMx}
      done

    elif [ ${run_num} == "3" ]; then
      echo "******************************************"
      echo "Running acf with 3 concatenated runs"
      echo "******************************************"
      if [ ${fold_id} == "a" ]; then
        echo "** Running Fold 1: 3a **"
        for run in run2 run3 run4 
        do
        run3dFWHMx="3dFWHMx -mask ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_${run}/desc-optcom_bold_mask.nii.gz"
        run3dFWHMx="${run3dFWHMx} -input ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version_runs}/${subject}_${task}_${run}_errts.nii.gz -acf ${output_dir}/${subject}_${task}_${run}_acf_emp.1D"
        run3dFWHMx="${run3dFWHMx} > ${output_dir}/${subject}_${task}_${run}_acf.1D"
        eval ${run3dFWHMx}
        done
      elif [ ${fold_id} == "b" ]; then
        echo "** Running Fold 2: 3b **"
        for run in run1 run3 run4 
        do
        run3dFWHMx="3dFWHMx -mask ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_${run}/desc-optcom_bold_mask.nii.gz"
        run3dFWHMx="${run3dFWHMx} -input ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version_runs}/${subject}_${task}_${run}_errts.nii.gz -acf ${output_dir}/${subject}_${task}_${run}_acf_emp.1D"
        run3dFWHMx="${run3dFWHMx} > ${output_dir}/${subject}_${task}_${run}_acf.1D" 
        eval ${run3dFWHMx}
        done
      elif [ ${fold_id} == "c" ]; then
        echo "** Running Fold 3: 3c **"
        for run in run1 run2 run4 
        do
        run3dFWHMx="3dFWHMx -mask ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_${run}/desc-optcom_bold_mask.nii.gz"
        run3dFWHMx="${run3dFWHMx} -input ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version_runs}/${subject}_${task}_${run}_errts.nii.gz -acf ${output_dir}/${subject}_${task}_${run}_acf_emp.1D"
        run3dFWHMx="${run3dFWHMx} > ${output_dir}/${subject}_${task}_${run}_acf.1D" 
        eval ${run3dFWHMx}
        done
      elif [ ${fold_id} == "d" ]; then
        echo "** Running Fold 4: 3d **"
        for run in run1 run2 run3
        do
        run3dFWHMx="3dFWHMx -mask ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.tedana_${run}/desc-optcom_bold_mask.nii.gz"
        run3dFWHMx="${run3dFWHMx} -input ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version_runs}/${subject}_${task}_${run}_errts.nii.gz -acf ${output_dir}/${subject}_${task}_${run}_acf_emp.1D"
        run3dFWHMx="${run3dFWHMx} > ${output_dir}/${subject}_${task}_${run}_acf.1D" 
        eval ${run3dFWHMx}
        done
      fi
    fi

    # Calculate average acf across all subjects/scans
    awk 'FNR == 1 { nfiles++; ncols = NF }
      { for (i = 1; i < NF; i++) sum[FNR,i] += $i
        if (FNR > maxnr) maxnr = FNR
      }
      END {
          for (line = 2; line <= maxnr; line++)
          {
              for (col = 1; col < ncols; col++)
                    printf "  %f", sum[line,col]/nfiles;
              printf "\n"
          }
      }' ${output_dir}/sub-*_acf.1D > ${output_dir}/allSub_acf.1D

  fi

  ##########################
  ### Cluster simulation ###
  ##########################

  if [ "${DO_ClustSim}" -eq 1 ]
  then
    echo "Cluster Simulation"
    echo "*****************"

    # Obtain ACF values
    acf_file="${output_dir}/allSub_acf.1D"

    read -ra values < ${acf_file}

    # Assign the values to individual variables
    value1=${values[0]}
    value2=${values[1]}
    value3=${values[2]}

    # Get cluster information
    3dClustSim -prefix "${output_dir}/${fxn}_${cond1}-${cond2}_${version}_clustSim" \
      -mask ${brain_mask} -acf  ${value1}  ${value2}  ${value3} -iter 10000

  fi

  ##########################
  ####### Clusterize #######
  ##########################

  if [ "${DO_Clusterize}" -eq 1 ]
  then
    echo "Clusterize"
    echo "*****************"

    # Find number of voxels at set pthr and athr =.05 (extract from clustSim file)
    clustSim_file="${output_dir}/${fxn}_${cond1}-${cond2}_${version}_clustSim.NN1_1sided.1D"
    target_column=3 #a = 0.05
    clustNum=$(awk -v row="$target_row" -v col="$target_column" '/^[^#]/{line_count++; if(line_count==row){for(i=1;i<=NF;i++) if($i!~/^#/){if(++count==col){print $i;exit}}}}' "$clustSim_file" | sed 's/^#.*$//')
    rounded_clustNum=$(( (${clustNum%.*} + 1) ))
    echo "Minimum cluster size is $rounded_clustNum"
    echo "Thresholded at p<0.${pthr} and clustered using alpha<0.${athr}"

    ## Threshold and cluster group results based on info from 3dClustSim
    ## Manually change the parameters below based on results of 3dClustSim and the desired clustering (e.g. -1sided, -2sided, -bisided)
    if [ ${run_num} == "4" ]; then
      echo "******************************************"
      echo "Running clusterize with 4 concatenated runs"
      echo "******************************************"
        
        input_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_ConsOrth_concat4runs_sm_1ses/${subject}_${task}_TaskReg_tstat.nii.gz"
    
    elif [ ${run_num} == "3" ]; then
      echo "******************************************"
      echo "Running clusterize with 3 concatenated runs"
      echo "******************************************"
        input_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_ConsOrth_concat3${fold_id}runs_sm_1ses/${subject}_${task}_TaskReg_tstat.nii.gz"
    
    fi

    3dClusterize -inset ${input_file} \
    -mask ${brain_mask} -ithr 0 -idat 0 -1sided RIGHT_TAIL p=0.${pthr} -NN 1 -clust_nvox ${rounded_clustNum} \
    -pref_dat "${output_dir}/${subject}_${fxn}_${cond1}-${cond2}_${model}_${version}_clusters_bcoef_p${pthr}_a${athr}.nii.gz"

  fi

  ##########################
  ### Func2Anat ###
  ##########################

  if [ "${DO_func2anat}" -eq 1 ]
  then
    echo "Func2Anat"
    echo "*****************"

    output_file="${output_dir}/${subject}_${fxn}_${cond1}-${cond2}_${model}_${version}_clusters_bcoef_p${pthr}_a${athr}"
      if [ ! -f ${output_file}_func2anat.nii.gz ] ; then
            # Transform output file (topPer) to anatomical space
            ./x.PreProc_Transform_lin.sh ${output_file} \
              ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
              ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
              ${output_dir} \
              func2anat

            fslmaths ${output_file}_func2anat.nii.gz -bin ${subject}_Cluster_${cond1}-${cond2}_clusters_bcoef_p${pthr}_a${athr}_func2anat_bin.nii.gz # binarize mask
      fi
    
    tstat_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_tstat"
      if [ ! -f ${tstat_file}_func2anat.nii.gz ] ; then
            # Transform output file (topPer) to anatomical space
            ./x.PreProc_Transform_lin.sh ${tstat_file} \
              ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
              ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
              ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
              func2anat
      fi

    bcoef_file="${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version}/${subject}_AUDIO_TaskReg_bcoef"
      if [ ! -f ${bcoef_file}_func2anat.nii.gz ] ; then
            # Transform output file (topPer) to anatomical space
            ./x.PreProc_Transform_lin.sh ${bcoef_file} \
              ${parent_dir}/derivatives/${subject}/ses-01/anat/T1_biascorr_brain \
              ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.reg_2mm/${subject}_${session}_${task}_func2anat \
              ${parent_dir}/derivatives/${subject}/${session}/func/${task}/output.GLM_${model}_${version} \
              func2anat
      fi
    
  fi

done ## end condition loop
done ## end subject loop
done ## end p-threshold loop