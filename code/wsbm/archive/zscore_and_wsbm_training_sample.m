%Running on the cluster
datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation'
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
z_avg_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'
addpath(genpath('/data/picsl/mackey_group/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'))
addpath(genpath('/home/utooley/matlab/WSBM_v1.2'))
addpath(genpath('/home/utooley/matlab/BCT'))
addpath(genpath('/home/utooley/matlab/NCT'))
addpath(genpath('/home/utooley/matlab/'))

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400Networks'
z_avg_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'

%get the subject list,excluding those who have NAs
%subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id(n)) %look at this
    %sub=subjlist{n,:}; 
    run1=char(subjlist.var1(n)) %%% FOR RUN 1 %%
    file=fullfile(datadir,sub, run1, strcat('fcon/schaefer400x7/',sub,'_',run1,'_schaefer400x7_network.txt'));
    try
    %subfcmat=load(file{1});
    subfcmat=load(file);
    %make into adjacency matrix and save out
    size_vec=tril(ones(400,400),-1);
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_',run1,'_schaefer400x7_network.txt'));
    csvwrite(outfile, subfcmat);
   % subfcmat(:,103)=[]; %never checked parcel coverage for this.
    %replace the diagonal of 1's with 0's
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc1=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc1(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(sub,'_', run1,'_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, zfc1);
    run2=char(subjlist.var2(n))  %%% FOR RUN 2 %%
    file=fullfile(datadir,sub, run2, strcat('fcon/schaefer400x7/',sub,'_',run2,'_schaefer400x7_network.txt'));
    subfcmat=load(file);
    size_vec=tril(ones(400,400),-1); %make into adjacency matrix and save out
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_',run2,'_schaefer400x7_network.txt'));
    csvwrite(outfile, subfcmat);
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc2=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc2(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(sub,'_', run2,'_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, zfc2);
    %% AVERAGE THE TWO RUNS %%
    aggregmat=(zfc2+zfc1)/2;
    outfile=fullfile(z_avg_outdir, strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, aggregmat);
    catch
    fprintf('Cant read sub %s run %s, skipped. \n', sub);
  end
end


%% Run WBSM model signed matrices with different values of k
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/search_over_k'
%on cluster
wsbm_dir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/search_over_k'
%add a loop for different values of k (communities)?
%for k=1:12, for example
%or when running on the cluster can use the function below
%wsbm_function(sub, comm) %need to set outputs in the WSBM function so that it writes to the right place.
%preallocate for speed

k=[4:10] %change this to be the right range
names=strcat('k',string(k))
for k=6:7 %change this to be the right range
    num_comms=k+3
    sub_log_lik_for_weights.(names{k})=zeros(height(subjlist),1);
    sub_log_lik_for_edges.(names{k})=zeros(height(subjlist),1);
    sub_log_evidence.(names{k})=zeros(height(subjlist),1);
    sub_labels.(names{k})=zeros(height(subjlist),400);
    temp_weights=zeros(height(subjlist),1);
    temp_edges=zeros(height(subjlist),1);
    temp_evidence=zeros(height(subjlist),1);
    temp_labels=zeros(height(subjlist),400);
    %parfor n=1:height(subjlist)
    parfor n=1:height(subjlist)
        sub=char(subjlist.id(n)) %look at this
        [temp_weights(n,1),temp_edges(n,1),temp_evidence(n,1),temp_labels(n,:)]=wsbm_function(sub, num_comms,wsbm_dir,z_avg_outdir)
    end
    sub_log_lik_for_weights.(names{k})=temp_weights
    sub_log_lik_for_edges.(names{k})=temp_edges
    sub_log_evidence.(names{k})=temp_evidence
    sub_labels.(names{k})=temp_labels
    
       outfile1=dataset(subjlist.id, sub_log_lik_for_weights.(names{k}), sub_log_lik_for_edges.(names{k}), sub_log_evidence.(names{k}), sub_labels.(names{k}))
       export(outfile1,'File',fullfile(wsbm_dir,strcat('wsbm_search_over_k', num2str(num_comms),'_n670_site16_30trials.csv')),'Delimiter',',')
    end

%save subject-level variables.
outfile=dataset(subjlist, sub_log_lik_for_weights, sub_log_lik_for_edges, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_search_over_k_n670_site16_30trials.csv'),'Delimiter',',')

%% Relabel community partitions to be most parsimonious across subjects, create a consensus partition
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
yeo_nodes=dlmread('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

%read in each subject's wsbm partition
for n=1:height(subjlist)
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'));
    try 
        load(file);
        sub_log_evidence(n,1)=Model.Para.LogEvidence;
        %a row for each subject's labels for the k=7 partition
        %sub_labels(n,:)=Labels; skip this, before relabeling for
        %persistence with Yeo these are not right
        %Look at how many unique labels were output, since Rick says these might be
        %different
        sub_unique_labels(n,1)=numel(unique(Labels));
        models(n, :)=Model;
        partition=Labels;%load in the partition from the WSBM for this subject
        %relabel WSBM to fit numbers label to Yeo partition so can compare stats
        %relabel for similarity to Yeo partition
        %temp=multislice_pair_labeling([yeo_nodes partition]);
        %sub_labels(n,:)=temp(:,2)';
        %partition=temp(:,2);
        %generate one big partition matrix
        part_matrix(:,n)=partition;
        %save how similar a subject is to Yeo partition
        sub_similarity_to_yeo=zrand(partition, yeo_nodes);
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%save subject-level variables.
outfile=dataset(subjlist.id, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_n670_site16_50trials.csv'),'Delimiter',',')

%either use multislice_pair_label or just pair_label
%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
%don't do this if you relabel for Yeo at the subject level!
yeo_opt_part_matrix=multislice_pair_labeling([yeo_nodes part_matrix])
noyeo_opt_part_matrix=multislice_pair_labeling(opt_part_matrix)
%opt_part_matrix=part_matrix

node_var=var(opt_part_matrix, 0,2) %calculate node-wise SD (how much the labeling varies across subjects after each is relabeled for Yeo)
%node_entropy=Entropy(opt_part_matrix') %first entropy from Entropy.m
h = hist(opt_part_matrix(:,:)',7); % node histogram
p = bsxfun(@rdivide,h,sum(h));  % probabilities
node_entropy = -nansum(p.*log2(p))';        % entropy

[node_mode freq ties]=mode(opt_part_matrix, 2)%calculate nodewise mode (most common label), and export ties as well.

%% Create a consensus partition
consensus_mat=consensus_similarity(opt_part_matrix')'; %maybe transposed?
%see z-score of the Rand coefficient for the consensus community
[zrandconsensus,~,~,VIconsensus]=zrand(yeo_nodes,consensus_mat)
[zrandmodal,~,~,VImodal]=zrand(yeo_nodes,node_mode)
%% Permutation testing for ZRand
for i=1:10000
    permuted_labels=consensus_mat(randperm(length(consensus_mat)));
    zrand_permuted(i)=zrand(yeo_nodes, permuted_labels');
end
hist(zrand_permuted)
hold on
line([zrandconsensus zrandconsensus], [0,300])
perm_zrand_dist=fitdist(zrand_permuted','Normal')
p=cdf(perm_zrand_dist, zrandconsensus)
%% SAVE OUTFILES
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/'
outfile=(fullfile(outdir, 'n670_training_sample_relabel_first_consensus_mat_and_nodal_variance.mat'))
save(outfile, 'opt_part_matrix', 'part_matrix', 'node_var', 'node_entropy', 'consensus_mat', 'node_mode')
outfile=(fullfile(outdir, 'n670_training_sample_relabel_first_nodal_var.mat'))
save(outfile, 'node_var')
outfile=(fullfile(outdir, 'n670_training_sample_relabel_first_nodal_entropy.mat'))
save(outfile, 'e')
outfile=(fullfile(outdir, 'consensus_mat.mat'))
save(outfile, 'consensus_mat')

%transpose the re-labeled partition matrix so it can be saved with each
%subj on a row
% sub_labels_opt=opt_part_matrix'
% outfile=dataset(subjlist, sub_log_evidence, sub_unique_labels, sub_labels_opt)

%% Save a .nii file so can see WSBM on the brain 
%read into nifti?
templateVolume = '/data/picsl/mackey_group/tools/schaefer400/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'
nii = load_nii(templateVolume);
image = double(nii.img);
spacing = nii.hdr.dime.pixdim(2:4);

%% Assign labels to brain regions
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/'
outfile=(fullfile(outdir, 'n670_training_sample_relabel_first_consensus_mat_and_nodal_variance.mat'))
load(outfile)
% read in the mapping of template nifti voxel label numbers to actual brain regions
mapping = readtable('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7NodeIndex.1D.txt','ReadVariableNames',false);
%mapping = readtable('~/Documents/tooleyEnviNetworks/parcels/Glasser/glasser_lookup.csv');
mapping_nums=table2array(mapping(:,1));
brains=struct('node_var',node_var, 'node_entropy', node_entropy, 'node_mode', node_mode, 'consensus_mat', consensus_mat);
names=fieldnames(brains);
for x=1:numel(names)
    parc = zeros(size(image)); %create a new image matrix
    labels=brains.(names{x})
    for r=1:length(labels) %go through all possible voxel values and assign labels
        parc(find(image==mapping_nums(r))) = labels(r);
    end
    % WRITE OUTFILE
    parc = double(parc);
    orig = nii.hdr.hist.originator; %get the origin of the original image
    orig = orig(1:3);
    niiNew = make_nii(parc,spacing,orig); %write out the new nifti
    niiNew.hrd.dime.bitpix=16; %set the datatype
    save_nii(niiNew,strcat('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/wsbm_k7_n670_site16_',char(names(x)),'.nii'));
end
%% Look at node variance & entropy by Yeo network
%node_var=load(fullfile(dir,'node_var.mat'));
comms=unique(yeo_nodes)
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = node_var(Wi,1); % extract node var for those parcels
    %average it
    var_by_yeocomm(i)=mean(Wv_temp); 
end
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = node_entropy(Wi,1); % extract node var for those parcels
    %average it
    entropy_by_yeocomm(i)=mean(Wv_temp); 
end
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/'
export(dataset(var_by_yeocomm, entropy_by_yeocomm),'File',fullfile(outdir,'variance_entropy_in_wsbm_training_sample_assign_by_yeonet.csv'),'Delimiter',',')


