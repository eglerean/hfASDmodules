clear all
close all

load scores % we have behav scores for AS, as well as IDs for AS and NTs

% gather all framewise displacements
basepath=['/triton/becs/scratch/braindata/eglerean/ABIDE/00'];
Nsubj=27*2;
all_mean_FD=zeros(Nsubj,1);
all_var_FD=zeros(Nsubj,1);
FDth=zeros(Nsubj,2); % timepoints under 0.2 TH, and all number of time points
TRs=[2 2 1.5]; %CALTECH, CMU, PITT
all_subj_ids=[NT_ids; AS_ids];
sites_ids=[sites_ids;sites_ids];

% for all NT
for s=1:Nsubj;
    subjects{s,1}=[basepath num2str(all_subj_ids(s))];
    subject_TR(s)=TRs(sites_ids(s));
	getrid=20;
	if(subject_TR(s)==2) getrid=10; end
	temp=load([subjects{s,1} '/bramila/cfg.mat']);
	tempFD=temp.cfg.fDisplacement;
	tempFD(1:getrid)=[];
	tempFD(end+1-getrid:end)=[];
	all_mean_FD(s,1)=mean(tempFD);
	all_var_FD(s,1)=var(tempFD);
	FDth(s,1)=length(find(tempFD<0.5));
	FDth(s,2)=length(tempFD);
end
