close all
clear all
addpath('cbrewer'); % Cbrewer colormaps http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab

load scores
load si_median
load fd_mean
subnetwork_labels={'DM','LAN','AUD','SAL','Parietal','DA','SM','V1','VTL','Precuneus','Cerebellum','VIS'};

NPERMS=1e6;	% Number of permutations
% output range of scores
for s=1:length(score_labels)
	fprintf([score_labels{s} ': ' num2str(min(scores(:,s))) ' - ' num2str(max(scores(:,s))) '; '])
end
disp(' ');
% make models based on ADI-R/ADOS scores
temp=squareform(pdist(zscore(scores),'euclidean'));
model=1-temp/max(temp(:)); % 27 x 27 matrix of similarity between ABIDE subject pairs based on ADI-R/ADOS scores

temp=squareform(pdist(fd_mean,'euclidean'));
fd_model=1-temp/max(temp(:)); % 27 x 27 matrix of similarity between ABIDE subject pairs based on average framewise displacement

% NOTE: for old version of Matlab, you need to manually start matlabpool to use the parallel code
[motion_corr motion_p]=bramila_mantel(model,fd_model,NPERMS,'pearson'); % results for head motion similarity versus ADI-R/ADOS scores similarity

% for each NT subnetwork
for mod=1:12
    si=si_median(:,:,mod); % load precomputed SI (done in python)
	[cm(mod) pm(mod)]=bramila_mantel(model,si,NPERMS,'pearson'); % compute mantel test between subnetwork SI and model simm matrix based on ADI-R/ADOS
	[cm_mot(mod) pm_mot(mod)]=bramila_mantel(fd_model,si,NPERMS,'pearson');
	
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
	axis([1 1 0 M])
	axis square
	xlabel('Similarity of ADI-R/ADOS for subject-pair (ABIDE)')
	ylabel(['Similarity of ' subnetwork_labels{mod} ' subnetwork for subject-pair (ABIDE)'])
	set(gcf,'color',[1 1 1])
	set(gcf,'units','normalized','outerposition',[0 0 1 1]*.5)
	set(gca,'FontName','Arial')
	set(gca,'FontSize',10)
	saveas(gcf,['ADIRADOS_vs_SI_' num2str(mod) '.eps'],'psc2')
end
save ABIDE_results cm pm model fd_model motion_corr motion_p cm_mot pm_mot



