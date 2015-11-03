clear all
close all

load('/triton/becs/scratch/braindata/eglerean/ABIDE/git/hfASDmodules/ABIDE/scores.mat');
load('/triton/becs/scratch/braindata/eglerean/ABIDE/git/hfASDmodules/ABIDE/fd_mean_NTAS.mat')
load('/triton/becs/scratch/braindata/eglerean/ASnets/git/hfasdmodules/motion/AS_blacklist.mat');
ids=[NT_ids; AS_ids];
for s=1:length(ids)
    data=bramila_loadpickle(['/triton/becs/scratch/braindata/eglerean/ABIDE/Brain.Asperger.Rest/Data_v3/louvain_clusters_00' num2str(ids(s)) '_density_0.02_with_mst.pkl']);
    % find best module
    best=find(max(data.modularity_q)==data.modularity_q);
    best=best(1); % in case two are the same
    clust=data.louvain_clusters;
    clust(:,blacklist)=[]; % get rid of blacklist
    if(any(isnan(clust(:)))) error('there shouldnt be nans'); end
    best_mod(s,:)=[data.modularity_q(best) length(unique(clust(best,:)))];
    nmods=max(clust+1,[],2);
    avg_mod(s,:)=[mean(data.modularity_q) std(data.modularity_q) mean(nmods) std(nmods)];
end

% check the best module modularity and number of modules versus mean framewise displacement, std of framewise displacement
[cm_best_modula pm_best_modula]=corr(fd_mean,best_mod(:,1),'type','Spearman')
[cm_best_modnum pm_best_modnum]=corr(fd_mean,best_mod(:,2),'type','Spearman')

% extra check: see if average modularity/module num is correlated to mean framewise displacement
[cm_avg pm_avg]=corr(fd_mean,avg_mod,'type','Spearman')

