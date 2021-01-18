#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/spatial_topography/output/job_output/
#$ -l h_vmem=7.5G,s_vmem=7.3G
#$ -V
#$ -cwd

#Note that you have to do the last subject by hand! The cvs does not read the last subject.

subjlist=/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n670_filtered_runs_site16_postprocess_for_snr.csv
#subjlist=/data/jux/mackey_group/Ursula/n2_filtered_runs_site16_postprocess.csv
while IFS=, read -r col1 sub run1 run2 col5 col6
do
    echo "I got:$col1|$sub|$run1|$run2"
    sub_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/${sub}
    #truncate run1
    fslroi ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed.nii.gz ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated.nii.gz 0 370
    #3dcalc -a ${sub_dir}/${sub}_task-rest_${run1}_space-T1w_desc-preproc_bold.nii.gz[0-369] -expr 'a' -prefix ${sub_dir}/${sub}_task-rest_${run1}_space-T1w_desc-preproc_bold_truncated.nii.gz
    #might be faster ^
    #truncate run2
    fslroi ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed.nii.gz ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated.nii.gz 0 370
    #3dcalc -a ${sub_dir}/${sub}_task-rest_${run2}_space-T1w_desc-preproc_bold.nii.gz[0-369] -expr 'a' -prefix ${sub_dir}/${sub}_task-rest_${run2}_space-T1w_desc-preproc_bold_truncated.nii.gz
    echo first run
    fslmaths ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated.nii.gz -Tmean ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_mean.nii.gz #input the fmriprep run, but first truncate it to 370 volumes.
    fslmaths ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated.nii.gz -Tstd ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_std.nii.gz
    fslmaths ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_mean.nii.gz -div ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_std.nii.gz ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_tSNR.nii.gz
    echo second run
    fslmaths ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated.nii.gz -Tmean ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_mean.nii.gz #input the fmriprep run, but first truncate it to 370 volumes.
    fslmaths ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated.nii.gz -Tstd ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_std.nii.gz
    fslmaths ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_mean.nii.gz -div ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_std.nii.gz ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_tSNR.nii.gz
    #average run 1 and run 2
    fslmaths ${sub_dir}/${run1}/prestats/${sub}_${run1}_preprocessed_truncated_tSNR.nii.gz -add ${sub_dir}/${run2}/prestats/${sub}_${run2}_preprocessed_truncated_tSNR.nii.gz -div 2 ${sub_dir}/${sub}_preprocessed_truncated_tSNR_2run_avg.nii.gz
done < $subjlist

# #Transform tsnr maps to fsaverage6 space for each participant, average across them
#
# export SUBJECTS_DIR=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/freesurfer
# data_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/fmriprep
# out_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/surfaces
#
# subject=${1}
# fs_subject=${subject}
#
# sub_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/fmriprep/${sub}/func/
#
# for sub in
# do
# if [ -f ${out_dir}/${subject}/surf/lh.fs5_${subject}_task-rest_tsnr_avg_2runs_trunc_afterfmriprep.nii.gz ]; then
#   echo $sub $sesh 'done already'
#   # if it exists already, skip it
# else
#
# for hem in lh rh
# do
#   echo hem is $hem
# sub_dir=/cbica/projects/spatial_topography/dropbox/data/${sub}/func/
# SUBJECTS_DIR=/cbica/projects/spatial_topography/dropbox/data/freesurfer/freesurfer
# mri_vol2surf --src ${sub_dir}/${sub}_task-rest_space-T1w_desc-preproc_bold_truncated_tSNR_2run_avg.nii.gz \
#               --out ${hem}.${sub}_task-rest_space-T1w_desc-preproc_bold_truncated_tSNR_2run_avg.nii.gz \
#               --regheader ${sub} \
#               --trgsubject fsaverage5 \ #NEED TO CHECK THAT FSAVERAGE5 LOOKS NICE STILL!
#               --hemi ${hem}
# done

#this looks good, check this by running freeview --recon ${sub} and overlaying the original functional volume on the T1.mgz, and the functional image on the lh.inflated of the subject.


#Then average across them.
