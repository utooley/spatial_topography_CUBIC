

#LH
#upsample the fsaverage5 confidence maps to fsaverage6 space
annot_file_in=/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/lh.silhouette.fsaverage5.curv
annot_file_out=/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/lh.silhouette.fsaverage6.curv

mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage6 --hemi lh --sval $annot_file_in --src_type curv --tval $annot_file_out --trg_type curv

#RH
annot_file_in=/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/label/rh.silhouette.fsaverage5.curv
annot_file_out=/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage6/label/rh.silhouette.fsaverage6.curv

mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage6 --hemi rh --sval $annot_file_in --src_type curv --tval $annot_file_out --trg_type curv

#LH
#upsample the fsaverage5 snr maps to fsaverage6 space
snr_file_in=~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/lh.1000subjects.snr.nii.gz
snr_file_out=~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/lh.fs6.1000subjects.snr.nii.gz
mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage6 --hemi lh --sval $snr_file_in --tval $snr_file_out --trg_type curv
#RH
snr_file_in=~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/rh.1000subjects.snr.nii.gz
snr_file_out=~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/rh.fs6.1000subjects.snr.nii.gz

mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage6 --hemi rh --sval $snr_file_in --tval $snr_file_out --trg_type curv
