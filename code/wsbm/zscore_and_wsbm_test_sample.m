% %Running on the cluster
% datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
% listdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation'
% outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
% z_avg_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
% wsbm_dir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'
% addpath(genpath('/data/picsl/mackey_group/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'))
% addpath(genpath('/home/utooley/matlab/WSBM_v1.2'))
% addpath(genpath('/home/utooley/matlab/BCT'))
% addpath(genpath('/home/utooley/matlab/NCT'))
% addpath(genpath('/home/utooley/matlab/'))

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'
z_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zNetworks'
noz_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400Networks'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample'
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/NCT'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/'))


%get the subject list,excluding those who have NAs
%subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=readtable(fullfile(listdir,'n544_filtered_runs_site14site20_postprocess.csv')) %the first two runs here are those that were input into gwMRF

%% Z-score FC matrices
% for n=1:height(subjlist)
%     sub=char(subjlist.id(n)) %look at this
%     %sub=subjlist{n,:}; 
%     run1=char(subjlist.var1(n)) %%% FOR RUN 1 %%
%     file=fullfile(datadir,sub, run1, strcat('fcon/schaefer400x7/',sub,'_',run1,'_schaefer400x7_network.txt'));
%     try
%     %subfcmat=load(file{1});
%     subfcmat=load(file);
%     %make into adjacency matrix and save out
%     size_vec=tril(ones(400,400),-1);
%     adj_mat=size_vec;
%     adj_mat(adj_mat==1)=subfcmat;
%     subfcmat=adj_mat+adj_mat';
%     outfile=fullfile(noz_outdir,strcat(sub,'_',run1,'_schaefer400x7_network.txt'));
%     csvwrite(outfile, subfcmat);
%    % subfcmat(:,103)=[]; %never checked parcel coverage for this.
%     %replace the diagonal of 1's with 0's
%     for x=1:400
%         subfcmat(x,x)=0;
%     end
%     %create an empty z-matrx
%     zfc1=[];
%     for i=1:400
%         %cycle through each column of the FC matrix and do a fisher r-to-z
%         %for each value
%         zfc1(:,i)=fisherz(subfcmat(:,i));
%     end
%     outfile=fullfile(z_outdir, strcat(sub,'_', run1,'_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, zfc1);
%     run2=char(subjlist.var2(n))  %%% FOR RUN 2 %%
%     file=fullfile(datadir,sub, run2, strcat('fcon/schaefer400x7/',sub,'_',run2,'_schaefer400x7_network.txt'));
%     subfcmat=load(file);
%     size_vec=tril(ones(400,400),-1); %make into adjacency matrix and save out
%     adj_mat=size_vec;
%     adj_mat(adj_mat==1)=subfcmat;
%     subfcmat=adj_mat+adj_mat';
%     outfile=fullfile(noz_outdir,strcat(sub,'_',run2,'_schaefer400x7_network.txt'));
%     csvwrite(outfile, subfcmat);
%     for x=1:400
%         subfcmat(x,x)=0;
%     end
%     %create an empty z-matrx
%     zfc2=[];
%     for i=1:400
%         %cycle through each column of the FC matrix and do a fisher r-to-z
%         %for each value
%         zfc2(:,i)=fisherz(subfcmat(:,i));
%     end
%     outfile=fullfile(z_outdir, strcat(sub,'_', run2,'_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, zfc2);
%     %% AVERAGE THE TWO RUNS %%
%     aggregmat=(zfc2+zfc1)/2;
%     outfile=fullfile(z_avg_outdir, strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, aggregmat);
%     catch
%     fprintf('Cant read sub %s run %s, skipped. \n', sub);
%   end
% end
% 

%% Run WSBM replication, 100% density and k=7
k=7
%densities={'0.5','0.7'}
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'

% for m=1:length(densities)
%     density=densities{m}
for n=1:height(subjlist)
    sub=char(subjlist.id(n)) %look at this
    outfile=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'))
    if (exist(outfile)==2) %if it's already written don't do it again
        fprintf('Sub %s already exists. \n', sub);
    else
    file=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
    try
    subfcmat=load(file);
    [Labels Model]=wsbm(subfcmat, k,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 50, 'verbosity', 1, 'alpha', 0, 'parallel', 0);
    %logHw - Additive Log-likelihood constant for W_Distr
    sub_log_lik_for_weights(n,1)=Model.Data.logHw;
    %logHe - Additive Log-likelihood constant for E_Distr
    sub_log_lik_for_edges(n,1)=Model.Data.logHe;
    %LogEvidence - Marginal Log-likelihood (aka Log-Evidence), a model selection criterion 
    sub_log_evidence(n,1)=Model.Para.LogEvidence;
    %a row for each subject's labels for the k=7 partition
    sub_labels(n,:)=Labels;
    %Look at how many unique labels were output, since Rick says these might be
    %different
    sub_unique_labels(n,1)=numel(unique(Labels));
    models(n, :)=Model;
    %save the model
    save(outfile, 'Model', 'Labels')
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
    end
end
outfile=dataset(subjlist, sub_log_lik_for_weights, sub_log_lik_for_edges, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_n544_site14site20_50trials.csv'),'Delimiter',',')

%% Run WSBM with threshold at 50% and 70% density
% k=7
% densities={'0.5','0.7'}
% wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/'
% z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
% 
% for m=1:length(densities)
%     density=densities{m}
% for n=1:height(subjlist)
%     sub=char(subjlist.id(n)) %look at this
%     outfile=fullfile(wsbm_dir,strcat('thresh_',density,'/',sub,'_wsbm_thresh.mat'))
%     if (exist(outfile)==2) %if it's already written don't do it again
%         fprintf('Sub %s already exists. \n', sub);
%     else
%     file=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
%     try
%     subfcmat=load(file);
%     subfcmat=threshold_proportional(subfcmat, density);
%     [Labels Model]=wsbm(subfcmat, k,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 50, 'verbosity', 1, 'alpha', 0, 'parallel', 0);
%     %logHw - Additive Log-likelihood constant for W_Distr
%     sub_log_lik_for_weights(n,1)=Model.Data.logHw;
%     %logHe - Additive Log-likelihood constant for E_Distr
%     sub_log_lik_for_edges(n,1)=Model.Data.logHe;
%     %LogEvidence - Marginal Log-likelihood (aka Log-Evidence), a model selection criterion 
%     sub_log_evidence(n,1)=Model.Para.LogEvidence;
%     %a row for each subject's labels for the k=7 partition
%     sub_labels(n,:)=Labels;
%     %Look at how many unique labels were output, since Rick says these might be
%     %different
%     sub_unique_labels(n,1)=numel(unique(Labels));
%     models(n, :)=Model;
%     %save the model
%     save(outfile, 'Model', 'Labels')
%     catch
%     fprintf('Cant read sub %s, skipped. \n', sub);
%     end
%     end
% end
% outfile=dataset(subjlist, sub_log_lik_for_weights, sub_log_lik_for_edges, sub_log_evidence, sub_unique_labels, sub_labels)
% export(outfile,'File',fullfile(wsbm_dir,strcat('thresh_',density),'wsbm_k7_n670_site16_thresh_50trials.csv'),'Delimiter',',')
% end
% 
%% Relabel community partitions to be most parsimonious across subjects, create a consensus partition
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
yeo_nodes=dlmread('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
yeo_nodes=dlmread('/cbica/projects/spatial_topography/tools/parcellations/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

%read in each subject's wsbm partition
for n=1:height(subjlist)
    sub=char(subjlist.id(n))
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
        part_matrix(:,n)=partition;
        %save how similar a subject is to Yeo partition
        %sub_similarity_to_yeo=zrand(partition, yeo_nodes);
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%save subject-level variables.
outfile=dataset(subjlist.id, sub_log_evidence, sub_unique_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_n544_site14site20_50trials.csv'),'Delimiter',',')

%either use multislice_pair_label or just pair_label
%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
%don't do this if you relabel for Yeo at the subject level!
%yeo_opt_part_matrix=multislice_pair_labeling([yeo_nodes part_matrix])
opt_part_matrix=part_matrix
noyeo_opt_part_matrix=multislice_pair_labeling(opt_part_matrix)


%% Create consensus partitions
%Similarity, similarity partition
consensus_mat=consensus_similarity(noyeo_opt_part_matrix')'; %maybe transposed?
%see z-score of the Rand coefficient for the consensus community
[zrandconsensus,~,~,VIconsensus]=zrand(yeo_nodes,consensus_mat)

%consensus iterative partition (does it matter if you relabel or not? Not more than it varies with iterating over mod max)
gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(noyeo_opt_part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %can also look at nodes that still aren't able to be assigned , which nodes are have the most variance across nodes in the coocurrence matrix.
%see z-score of the Rand coefficient for the consensus iterative partition
zrandconsensus_iter=zrand(yeo_nodes,consensus_iter_mode)

%% Variance in assignment
%looking at the frequencies of assignments in consensus iterative across Yeo communities
freq_continuous=abs(freq-670)' %get the absolute value of how many times less than 670 it was assigned to the same community
comms=unique(yeo_nodes)
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = freq_continuous(Wi,1); % extract node var for those parcels
    %average it
    num_ties_consensus_iter_by_yeocomm(i)=mean(Wv_temp); 
end

%need to figure out entropy on the co-occurence matrix--not sure this is
%right!
p = bsxfun(@rdivide,cooccurence_matrix,sum(cooccurence_matrix));  % probabilities
p = bsxfun(@rdivide,cooccurence_matrix,670);  % probabilities
p2=p;
p2(p==1)=NaN;
node_entropy = -nansum(p.*log2(p))';        % entropy

%save these outfiles
outfile=('/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/consensus_iter_freq.mat')
save(outfile, 'freq')
%% Relabel the consensus representative and group partitions so they are visually comparable to Yeo
%relabel the consensus partitions so that they are visually comparable to Yeo
consensus_iter_mode_yeorelabeled=multislice_pair_labeling([yeo_nodes consensus_iter_mode']);
consensus_iter_mode_yeorelabeled=consensus_iter_mode_yeorelabeled(:,2);
consensus_represent_yeorelabeled=multislice_pair_labeling([yeo_nodes consensus_mat_noyeo])
consensus_represent_yeorelabeled=consensus_represent_yeorelabeled(:,2)

outfile=('/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/n544_test_sample_consensus_partitions_yeorelabeled.mat')
save(outfile, 'freq','consensus_iter_mode_yeorelabeled')

%% Permutation testing for ZRand
% For similarity partition
for i=1:10000
    permuted_labels=consensus_mat(randperm(length(consensus_mat)));
    zrand_permuted(i)=zrand(yeo_nodes, permuted_labels');
end
hist(zrand_permuted)
hold on
line([zrandconsensus zrandconsensus], [0,300])
perm_zrand_dist=fitdist(zrand_permuted','Normal')
p=cdf(perm_zrand_dist, zrandconsensus)

% For consensus partition

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
% 
% %% UNUSED-Look at node variance & entropy by Yeo network
% %node_var=load(fullfile(dir,'node_var.mat'));
% comms=unique(yeo_nodes)
% for i = 1:length(comms) % for each of the yeo communities
%     Wi = yeo_nodes == comms(i); % find index for parcels with that community
%     Wv_temp = node_var(Wi,1); % extract node var for those parcels
%     %average it
%     var_by_yeocomm(i)=mean(Wv_temp); 
% end
% for i = 1:length(comms) % for each of the yeo communities
%     Wi = yeo_nodes == comms(i); % find index for parcels with that community
%     Wv_temp = node_entropy(Wi,1); % extract node var for those parcels
%     %average it
%     entropy_by_yeocomm(i)=mean(Wv_temp); 
% end
% outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/'
% export(dataset(var_by_yeocomm, entropy_by_yeocomm),'File',fullfile(outdir,'variance_entropy_in_wsbm_training_sample_assign_by_yeonet.csv'),'Delimiter',',')


