%Running on the cluster
datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/'
z_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zNetworks'
noz_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400Networks'
z_avg_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'
addpath(genpath('/home/utooley/matlab/WSBM_v1.2'))
addpath(genpath('/home/utooley/matlab/system_matrix_tools/'))
addpath(genpath('/home/utooley/matlab/BCT/'))

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400Networks'
z_avg_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n546_filtered_runs_site14site20_postprocess.csv')) %the first two runs here are those that were input into gwMRF
subjlist=subjlist.id;

%% Run some partition quality metrics on the WSBM consensus partition and estimate network statistics
datadir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/wsbm/'
 %read in the WSBM partition
load('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling/n670_training_sample_consensus_partitions_yeorelabeled.mat')
consensus_nodes=consensus_iter_mode_yeorelabeled;
%preallocate variables
deviation_edge_weights_consensus=zeros(size(subjlist,1),1);
modul_consensus=zeros(size(subjlist,1),1);

%load in the models for each subject
for n=1:size(subjlist,1)
    sub=subjlist{n,:}
    try 
        %% Load FC matrix
        fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        %parcel 52 is already gone
        subfcmat = load(fcfile);
        %average edge weight for the subject
        avgweight(n,1)=mean(subfcmat(subfcmat~=0));
        %% Apply partition and calculate statistics with it
        [S, W, B] = segregation(subfcmat,consensus_nodes);
        system_segreg_consensus(n,1)=S;
        mean_within_sys_consensus(n,1)=W;
        mean_between_sys_consensus(n,1)=B;
        %Connectivity between each of the 7 blocks
        Ci=consensus_nodes;
        nCi = unique(Ci);
        M=subfcmat;

        for i = 1:length(nCi) % loop through communities
            for j = 1:length(nCi)
               Wi = Ci == nCi(i); % find index for within communitiy edges
               Bi = Ci == nCi(j); % find index for between communitiy edges to specific community

               Wv_temp = M(Wi,Wi); % extract within communitiy edges
               Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
                
               %calculate standard deviation of edges within blocks here
               deviation_edge_weights_consensus(n,1)=deviation_edge_weights_consensus(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
               %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
               Bv = [Bv_temp(:)'];
               system_connectivity(i,j)=mean(Bv(Bv~=0));
               %system_between(i,1)=mean(Bv);
               %if i==j
               %else
            end
        end

        %transpose these (or something) so they can be saved out on a subject basis.
        system_connectivity;
        system_connectivity_all(:,:,n)=system_connectivity;  %put this all into a matrix for everyone so that we can see the average
        consensus_conn_vector = reshape(system_connectivity',[],1)';
        consensus_connectivity(n,:)=consensus_conn_vector;
        %finish the std deviation of edge weights calculation
        deviation_edge_weights_consensus(n,1)=sqrt(deviation_edge_weights_consensus(n,1)/7);
        %some net stats on the WSBM partition
        %Participation coefficient average with consensus partition
        [Ppos Pneg]=participation_coef_sign(subfcmat, consensus_nodes);
        %each row is a subject, with part coef for 359 nodes
        sub_partcoef_pos_consensus(n,:)=mean(Ppos);
        sub_partcoef_neg_consensus(n,:)=mean(Pneg);
        %estimate modularity with this a-priori partition?
         %% CHECK THAT THESE TWO Q's ARE THE SAME %%
        [M Q]=modul_only(subfcmat, [], double(consensus_nodes), 'negative_asym');
        %QFModul(double(yeo_nodes), subfcmat) %this gives the same answer as modul_only above when using on a pos-only matrix, but 
        %we want the negative asymmteric weighting, so ignore for now.
        % [M Q]=community_louvain(subfcmat, [], double(yeo_nodes), 'negative_asym'); %this optimizes modularity, which is not what we ant.
        modul_consensus(n,1)=Q;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
end
%save outfile
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site14site20_test_sample'
outfile=dataset(subjlist, avgweight, system_segreg_consensus, mean_within_sys_consensus, mean_between_sys_consensus, deviation_edge_weights_consensus, sub_partcoef_pos_consensus, sub_partcoef_neg_consensus, modul_consensus, consensus_connectivity)
export(outfile,'File',fullfile(outdir,'n546_test_sample_schaefer400_wsbm_consensus_iter_network_stats.csv'),'Delimiter',',')

%also save the mean system connectivity matrix from wsbm
mean_system_conn_mat=mean(system_connectivity_all,3)
% header={'sys1', 'sys2', 'sys3','sys4','sys5','sys6','sys7'}
% mean_system_conn_mat.Properties.VariableNames=header;
save(fullfile(outdir, 'n546_mean_wsbm_consensus_iter_connectivity.mat'), 'mean_system_conn_mat')

%% Calculate log-evidence of different k's of WSBM calculated on testing sample
%idea is that higher k's could just be overfitting data, test the consensus
%partition on the replication sample

%% First must create a consensus partition at each k
yeo_nodes=dlmread('/cbica/projects/cbpd_main_data/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/search_over_k/'
%get the testing subject list
test_subjlist=readtable(fullfile(listdir,'n544_filtered_runs_site14site20_postprocess.csv')) %the first two runs here are those that were input into gwMRF
test_subjlist=test_subjlist.id;
%get training subject list
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

k=8
%read in each subject's wsbm partition
for k=[2,8,10,15]
for n=1:height(subjlist)
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_k',num2str(k),'_wsbm.mat'));
    try 
        load(file);
        partition=Labels;%load in the partition from the WSBM for this subject
        part_matrix(:,n)=partition;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%either use multislice_pair_label or just pair_label
%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
%don't do this if you relabel for Yeo at the subject level!
%yeo_opt_part_matrix=multislice_pair_labeling([yeo_nodes part_matrix])
opt_part_matrix=part_matrix
noyeo_opt_part_matrix=multislice_pair_labeling(opt_part_matrix)

%% Create consensus partitions

%consensus iterative partition (does it matter if you relabel or not? Not more than it varies with iterating over mod max)
gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(noyeo_opt_part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %can also look at nodes that still aren't able to be assigned , which nodes are have the most variance across nodes in the coocurrence matrix.
%save out
consensus_iter_mode_yeorelabeled=multislice_pair_labeling([yeo_nodes consensus_iter_mode']);
consensus_iter_mode_yeorelabeled=consensus_iter_mode_yeorelabeled(:,2);
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/'
outfile=(fullfile(outdir, strcat('n670_k',num2str(k),'training_sample_consensus_partition.mat')))
save(outfile, 'consensus_iter_mode_yeorelabeled')

%% Calculate log-evidence of this partition on the testing sample
%load subjlist
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'

%load in the FC matrix for each subject, estimate the log-evidence for the
%consensus partition at that k
for n=1:size(test_subjlist,1)
    sub=test_subjlist{n,:}
    try 
        %% Load FC matrix
        fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        %parcel 52 is already gone
        subfcmat = load(fcfile);
        %for connectivity matrix A
         Data = A(:);
        w = 0.000001;
        [WDistr] = setup_distr('Normal',[mean(Data),var(Data)+mean(Data)^2,w],1);

        %perform the following line just to get a basic model structure, run only
        %for a single trial
        [S(n,:),M] = wsbm(subfcmat,k,'W_Distr','Normal','numTrials',1,'E_Distr','None');
        LE(n) = calc_logEvidence(M);

        %want to get the log-likelihood for some partition L fitting the network A

        LP(n) = calc_nullEvidence(M,consensus_iter_mode_yeorelabeled,subfcmat);

    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
end

%save outfile
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/'
outfile=dataset(test_subjlist, LE,LP)
export(outfile,'File',fullfile(outdir,strcat('n544_test_sample_',num2str(k),'.csv')),'Delimiter',',')

end
