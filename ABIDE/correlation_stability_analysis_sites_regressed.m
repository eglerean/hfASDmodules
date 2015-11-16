clear all
close all
addpath('cbrewer');
map=cbrewer('qual','Set1',9);
% loading model and similarity matrices for AS group from ABIDE dataset
load scores
load ABIDE_results
load si_median

% result for the VTL subnetwork (the only one that is correlated with the model based on ADI-R/ADOS
vtl=si_median(:,:,9);
vtl_c=cm(9);
vtl_p=pm(9);


sites_model=1-squareform(pdist([sites_ids;sites_ids]));
sites_model(find(sites_model<0))=0;
sites_model_NT=sites_model(1:27,1:27);
sites_model_AS=sites_model_NT;

topids=find(triu(ones(size(vtl)),1));
    [tempA tempB tempC]=regress(vtl(topids),[sites_model(topids) ones(length(topids),1)]);
    vtl_sitesregressed=zeros(size(vtl));
    tempC=tempC-min(tempC); %normalizing similarity values
    vtl_sitesregressed(topids)=tempC/max(tempC);

    vtl_sitesregressed=vtl_sitesregressed+vtl_sitesregressed'+eye(size(vtl_sitesregressed));    % make it symmetrical
	vtl=vtl_sitesregressed;
	[vtl_c vtl_p]=bramila_mantel(vtl,model,1e4,'pearson');


% generate subsets of 13 to 25 subjects to run replication at multiple
% subsamples

rand('twister',0); % initialize the random generator
TotalNsubj=27;
NITER=351;  % 100 combinations for a subset of Nsubj subjects with Nsubj < TotalNsubj, 300 because "27 choose 25" = 351. 
for Nsubj=13:26;
	disp(['Considering a random subgroup of ' num2str(Nsubj) ' subjects']); 
   if(Nsubj==26)
        NITER=27; % if Nsubj is 26, then we only have (27 choose 26) = 27 combinations
    end
    seqs=zeros(NITER,Nsubj);    % number of combinations
    i=1;
    hash=-1*ones(NITER,1);      % hash table to check that the combination has not appeared before
    % while loop to generate the combinations
    while i<=NITER
        temp=randperm(TotalNsubj);
        temp=sort(temp(1:Nsubj));
        vec=zeros(TotalNsubj,1);
        vec(temp)=1;
        hash(i)=(2.^[1:27])*vec;
        if(length(find((hash(1:i-1)-hash(i))==0))>0)
            disp('this combination already exists');
        else
            seqs(i,:)=temp;
            i=i+1;
        end
    end


    for n =1:NITER; % for each combination compute mantel test correlation and p-value
		fprintf('%s',[num2str(round(10000*n/NITER)/100) '% '])
        subsa=seqs(n,:);
        modeltemp=model(subsa,subsa);
        vtltemp=vtl(subsa,subsa);
        [c(n) p(n)]=bramila_mantel(modeltemp,vtltemp,1e4,'pearson');
    end
	disp(' ')

    if(0) % disabled, only for testing purposes
        % plot for each Nsubj
        figure(Nsubj)
        subplot(2,1,1)
        hist(c,25);
        axis([0 0.45 0 25])
        hold on
        mm=median(c);
        plot([mm mm],[0 16],'r','LineWidth',2)
        subplot(2,1,2)
        semilogy(c,p,'k.')
        hold on
        plot([0 0.45],[ 0.05 0.05],'r')
        axis([0 0.45 0 1])
    end

    % final plot to see how correlation stabilizes with increasing Nsubj
    % For details see http://www.sciencedirect.com/science/article/pii/S0092656613000858
    figure(1000)
    hold on
    plot(Nsubj+.05*round(randn(size(c(find(p>0.05))))/.75),c(find(p>0.05)),'.','Color',map(9,:))
    plot(Nsubj+.05*round(randn(size(c(find(p<=0.05))))/.75),c(find(p<=0.05)),'.','Color',map(2,:))
    plot(Nsubj+[-.2 .2],[median(c) median(c)],'LineWidth',4,'Color',map(1,:))
    ppp(:,Nsubj-12)=prctile(c,[5 50 95]);
    
    % housekeeping
    clear c
    clear p
    clear hash
    clear seqs
end

figure(1000)
cosp=tanh(atanh(vtl_c)+0.1);
cosm=tanh(atanh(vtl_c)-0.1);
plot([12 27],[cosp cosp],'k--')
plot([12 27],[cosm cosm],'k--')
plot([12 27],[vtl_c vtl_c],'k')
plot(13:26,ppp(1,:),'k-.')
plot(13:26,ppp(3,:),'k-.')
axis([12 27 -.1 .4])
axis square
xlabel('Number of subjects')
ylabel('Correlation [Mantel Test]');
set(gcf,'Color',[1 1 1])
set(gcf,'units','normalized','outerposition',[0 0 1 1]*.5)
set(gca,'FontName','Arial')
set(gca,'FontSize',10)
saveas(gcf,['ABIDE_VTL_stability_sites_regressed.eps'],'psc2')
