close all
clear all
%% version 2 of script abide_mantel.m
% changes
%	- added 27 controls from ABIDE dataset
%	- check mantel test for mean framewise displacement for NT alone
%	- same as above for NT and AS combined
%	- check similarity based on site location for all, NT only, AS only
 

addpath('cbrewer'); % Cbrewer colormaps http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab
addpath('/triton/becs/scratch/braindata/shared/toolboxes/EffectSizeToolbox_v1.3/'); % http://www.mathworks.com/matlabcentral/fileexchange/32398-measures-of-effect-size-toolbox

load scores
load si_median_NTAS
load fd_mean_NTAS
subnetwork_labels={'DM','LAN','AUD','SAL','Parietal','DA','SM','V1','VTL','Precuneus','Cerebellum','VIS'};

NPERMS=1e6;	% Number of permutations
% output range of scores
for s=1:length(score_labels)
	fprintf([score_labels{s} ': ' num2str(min(scores(:,s))) ' - ' num2str(max(scores(:,s))) '; '])
end
disp(' ');

%% head motion models
% all subjects
temp=squareform(pdist(fd_mean,'euclidean'));
fd_model=1-temp/max(temp(:)); % 54 x 54 matrix of similarity between ABIDE subject pairs based on average framewise displacement
% NT only
temp=squareform(pdist(fd_mean(1:27),'euclidean'));
fd_model_NT=1-temp/max(temp(:));
% AS only
temp=squareform(pdist(fd_mean(28:end),'euclidean'));
fd_model_AS=1-temp/max(temp(:));

%% scanning sites models
sites_model=1-squareform(pdist([sites_ids;sites_ids]));
sites_model(find(sites_model<0))=0;
sites_model_NT=sites_model(1:27,1:27);
sites_model_AS=sites_model_NT;

%% model for permutation-based group comparisons
G=[ones(27) zeros(27); zeros(27) 2*ones(27)];
ids=find(triu(ones(54),1));
ttestmodel=G(ids);
goodids=find(ttestmodel>0);


for mod=1:12
	disp(num2str(mod)); 
	si=si_median(:,:,mod); % load precomputed SI (done in python)
	
	data=si(ids);
	stats=bramila_ttest2_np(data(goodids)',ttestmodel(goodids)',NPERMS);
	group_t(mod)=stats.tvals;
	group_p(mod)=min(stats.pvals);
	eff_si{mod}=mes(data(find(ttestmodel(goodids)==1)),data(find(ttestmodel(goodids)==2)),'hedgesg','nBoot',10000,'isDep',1);

	%% with sites regressed
	topids=find(triu(ones(size(si)),1));
    [tempA tempB tempC]=regress(si(topids),[sites_model(topids) ones(length(topids),1)]);
    si_sitesregressed=zeros(size(si));
    tempC=tempC-min(tempC); %normalizing similarity values
    si_sitesregressed(topids)=tempC/max(tempC);

    si_sitesregressed=si_sitesregressed+si_sitesregressed'+eye(size(si_sitesregressed));    % make it symmetrical
	data_sr=si_sitesregressed(ids);
	stats_sr=bramila_ttest2_np(data_sr(goodids)',ttestmodel(goodids)',NPERMS);
    group_t_sr(mod)=stats_sr.tvals;
    group_p_sr(mod)=min(stats_sr.pvals);
    eff_si_sr{mod}=mes(data_sr(find(ttestmodel(goodids)==1)),data_sr(find(ttestmodel(goodids)==2)),'hedgesg','nBoot',10000,'isDep',1);

	%% with sites regressed for two groups separately <- this is the one used in the paper
	si_NT=si(1:27,1:27);
	si_AS=si(28:end,28:end);
	topids=find(triu(ones(size(si_NT)),1));
    [tempA tempB tempC]=regress(si_NT(topids),[sites_model_NT(topids) ones(length(topids),1)]);
    si_sitesregressed_NT=zeros(size(si_NT));
    tempC=tempC-min(tempC); %normalizing similarity values
    si_sitesregressed_NT(topids)=tempC/max(tempC);
    si_sitesregressed_NT=si_sitesregressed_NT+si_sitesregressed_NT'+eye(size(si_sitesregressed_NT));    % make it symmetrical

	topids=find(triu(ones(size(si_AS)),1));
    [tempA tempB tempC]=regress(si_AS(topids),[sites_model_AS(topids) ones(length(topids),1)]);
    si_sitesregressed_AS=zeros(size(si_AS));
    tempC=tempC-min(tempC); %normalizing similarity values
    si_sitesregressed_AS(topids)=tempC/max(tempC);
    si_sitesregressed_AS=si_sitesregressed_AS+si_sitesregressed_AS'+eye(size(si_sitesregressed_AS));    % make it symmetrical

	% merging for within group, between group values are not used
	si_sitesregressed_wg=[si_sitesregressed_NT zeros(size(si_sitesregressed_NT)); zeros(size(si_sitesregressed_NT)) si_sitesregressed_AS];
	data_sr_wg=si_sitesregressed_wg(ids);
    stats_sr_wg=bramila_ttest2_np(data_sr_wg(goodids)',ttestmodel(goodids)',NPERMS);
    group_t_sr_wg(mod)=stats_sr_wg.tvals;
    group_p_sr_wg(mod)=min(stats_sr_wg.pvals);
    eff_si_sr_wg{mod}=mes(data_sr_wg(find(ttestmodel(goodids)==1)),data_sr_wg(find(ttestmodel(goodids)==2)),'hedgesg','nBoot',10000,'isDep',1);
	



	
	[cm_mot(mod) pm_mot(mod)]=bramila_mantel(fd_model,si,NPERMS,'pearson');
	[cm_mot_NT(mod) pm_mot_NT(mod)]=bramila_mantel(fd_model_NT,si_NT,NPERMS,'pearson');
	[cm_mot_AS(mod) pm_mot_AS(mod)]=bramila_mantel(fd_model_AS,si_AS,NPERMS,'pearson');

	[cm_site(mod) pm_site(mod)]=bramila_mantel(sites_model,si,NPERMS,'pearson');
	[cm_site_NT(mod) pm_site_NT(mod)]=bramila_mantel(sites_model_NT,si_NT,NPERMS,'pearson');
	[cm_site_AS(mod) pm_site_AS(mod)]=bramila_mantel(sites_model_AS,si_AS,NPERMS,'pearson');

end

% displays the results from the statistical tests and computes effect size
% with confidence interfals
disp(['Subnetwork ID|T-value|P-value|Hedges? g (95% c.i.)']);
for n=1:12;
    plot([n n],[eff_si_sr_wg{n}.hedgesgCi(1) eff_si_sr_wg{n}.hedgesgCi(2) ]);hold on;plot(n,eff_si_sr_wg{n}.hedgesg,'*');
    plot(n,-log10(group_p_sr_wg(n)),'r.');
    plot(n,group_t_sr_wg(n),'k*');
    stars='';
    if(group_p_sr_wg(n)<0.05/12) stars='*';end
    if(eff_si_sr_wg{n}.hedgesg>0.5) stars=[stars '*'];end
    disp([num2str(n) '|' num2str(group_t_sr_wg(n),'%1.3f') '|' num2str(group_p_sr_wg(n),'%0.4g') stars '|' num2str(eff_si_sr_wg{n}.hedgesg,'%1.3f') ' (' num2str(eff_si_sr_wg{n}.hedgesgCi(1),3) ' ' num2str(eff_si_sr_wg{n}.hedgesgCi(2),3) ')'])
end


% model based on ADI-R-ADOS scores (for AS subjects only)
temp=squareform(pdist(zscore(scores),'euclidean'));
model=1-temp/max(temp(:)); % 27 x 27 matrix of similarity between ABIDE subject pairs based on ADI-R/ADOS scores


% NOTE: for old version of Matlab, you need to manually start matlabpool to use the parallel code
[motion_corr motion_p]=bramila_mantel(model,fd_model_AS,NPERMS,'pearson'); % results for head motion similarity versus ADI-R/ADOS scores similarity
[sites_corr sites_p]=bramila_mantel(model,sites_model_AS,NPERMS,'pearson'); % results for head motion similarity versus ADI-R/ADOS scores similarity


% for each NT subnetwork
for mod=1:12
    disp(num2str(mod))
	si=si_median(:,:,mod); % load precomputed SI (done in python)
	si=si(28:end,28:end); % only AS
	
	[cm(mod) pm(mod)]=bramila_mantel(model,si,NPERMS,'pearson'); % compute mantel test between subnetwork SI and model simm matrix based on ADI-R/ADOS

	% regress out sites and build a new si matrix
	topids=find(triu(ones(size(si)),1));
	[tempA tempB tempC]=regress(si(topids),[sites_model_AS(topids) ones(length(topids),1)]);	
	si_sitesregressed=zeros(size(si));
	tempC=tempC-min(tempC);	%normalizing similarity values
	si_sitesregressed(topids)=tempC/max(tempC);
	
	si_sitesregressed=si_sitesregressed+si_sitesregressed'+eye(size(si_sitesregressed));	% make it symmetrical
	[cm_sitesregressed(mod) pm_sitesregressed(mod)]=bramila_mantel(model,si_sitesregressed,NPERMS,'pearson'); % compute mantel test between subnetwork SI and model simm matrix based on ADI-R/ADOS

	
	% plot mantel results
	figure(mod)
	map=cbrewer('qual','Set1',3);
	ids=find(triu(ones(size(si)),1)); % the top triangle of the similarity matrices for plotting
	plot(model(ids),si(ids),'.','Color',map(1,:),'MarkerSize',16);
	hold on
	ppp = polyfit((model(ids)), (si(ids)), 1);
	plot([0 1],polyval(ppp,[0 1]),'Color',map(1,:),'LineWidth',2)
	h=text(1,.05+polyval(ppp,1),['r = ' num2str(cm(mod),3) ]);
	set(h,'FontName','Arial')
	set(h,'HorizontalAlignment','right')
	h=text(1,.025+polyval(ppp,1),[' p = ' num2str(pm(mod),3) ]);
	set(h,'FontName','Arial')
	set(h,'HorizontalAlignment','right')
	uuu=unique(si(:));
	M=round(10*uuu(end-1))/10;
	axis([0 1 0 M])
	axis square
	xlabel('Similarity of ADI-R/ADOS for subject-pair (ABIDE)')
	ylabel(['Similarity of ' subnetwork_labels{mod} ' subnetwork for subject-pair (ABIDE)'])
	set(gcf,'color',[1 1 1])
	set(gcf,'units','normalized','outerposition',[0 0 1 1]*.5)
	set(gca,'FontName','Arial')
	set(gca,'FontSize',10)
	saveas(gcf,['ADIRADOS_vs_SI_' num2str(mod) '.eps'],'psc2')

	% plot mantel results for sites regressed
	figure(mod+100)
    map=cbrewer('qual','Set1',3);
    ids=find(triu(ones(size(si)),1)); % the top triangle of the similarity matrices for plotting
    plot(model(ids),si_sitesregressed(ids),'.','Color',map(1,:),'MarkerSize',16);
    hold on
    ppp = polyfit((model(ids)), (si_sitesregressed(ids)), 1);
    plot([0 1],polyval(ppp,[0 1]),'Color',map(1,:),'LineWidth',2)
    h=text(1,.05+polyval(ppp,1),['r = ' num2str(cm_sitesregressed(mod),3) ]);
    set(h,'FontName','Arial')
    set(h,'HorizontalAlignment','right')
    h=text(1,.025+polyval(ppp,1),[' p = ' num2str(pm_sitesregressed(mod),3) ]);
    set(h,'FontName','Arial')
    set(h,'HorizontalAlignment','right')
    uuu=unique(si_sitesregressed(:));
    M=round(10*uuu(end-1))/10;
    axis([0 1 0 M])
    axis square
    xlabel('Similarity of ADI-R/ADOS for subject-pair (ABIDE)')
    ylabel(['Similarity of ' subnetwork_labels{mod} ' subnetwork for subject-pair (ABIDE)'])
    set(gcf,'color',[1 1 1])
    set(gcf,'units','normalized','outerposition',[0 0 1 1]*.5)
    set(gca,'FontName','Arial')
    set(gca,'FontSize',10)
    saveas(gcf,['ADIRADOS_vs_SI_' num2str(mod) '_sitesregressed.eps'],'psc2')

end
save ABIDE_results_v2 cm pm model fd_model motion_corr motion_p cm_mot pm_mot



