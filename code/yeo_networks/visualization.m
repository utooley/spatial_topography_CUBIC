%Checking Yeo code outputs
cluster_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_2/'
file_name='yeo7_n670_2runsonly_1000tries_2.mat'

clustered = load(fullfile(cluster_dir,file_name));

ref_file = ['/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/1000subjects_clusters007_ref.mat'];
ref = load(ref_file); 

my_colors = [190   190  190
160    32  240
0     0  255
127   255  212
255   246  143
255     0  255
255   165    0
255     0    0]
    
    

%PLOT 7 NETWORKS-can change the colors to match WSBM colors?
CBIG_DrawSurfaceMaps(clustered.lh_labels, clustered.rh_labels, 'fsaverage6','inflated',0, 7, ref.colors)

CBIG_DrawSurfaceMaps(clustered.lh_labels, clustered.rh_labels, 'fsaverage5','inflated',0, 17, ref.colors)

%Silhouettes/confidence maps-can change the colors?
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/yeo7_n670_2runsonly_1000tries.mat')
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n131_5tries/silhou_cmap.mat')

CBIG_DrawSurfaceMaps(clustered.lh_s, clustered.rh_s, 'fsaverage6','pial',0, 0.4, color)

%read in WSBM data and relabel Schaefer400 annot file with WSBM assignments
%can use CBIG_read_annotation.m
%then use CBIG_VisualizeSurfaceAnnotationInFreeview to view.