% In this code, you will create fake 2D data and use clustering evaluation metrics
% to see how each will perform with respect to the fake data
% You can use real data in the form data = n x p where n is the observations and p are the
% variables.
% Fahd Alhazmi

% This block will prompt the user to create fake 2D points to evaluate them later on
 
clear all
clf
pos = get(gcf, 'position');
set(gcf, 'position',pos)
set(gca,'xlim',[0 1],'ylim',[0 1])
title(...
    {'Hit random positions to form clusters to test';...
    'press DELETE to remove last point';...
    'After adding all points, press ENTER'}...
    ,'FontSize',15)
[x, y] = getpts;
data=[x,y];
disp('Coordinates are now saved in ''artdata'' variable')



% create an output folder, unless there is already one

writepath= [pwd '/output'];
if ~exist(writepath)
    mkdir(writepath);
end

krange=[2:20];  % choose the range of clusters you want to evaluate
dist='sqEuclidean';  % choose the distance metric. See the help file for kmeans or other algorithms
runs=100;  % how many runs for each choice of clusters; the more the better CI are 
additionals=false;   % this will run additional similarity metrics at each k and will take additional time

[sMean,sStd,qualVal,qualStd,...
    subj_stab_all,...
    jc,jd,ri,vi,nmi]= ...
    eval_clustrs(data, krange, dist, runs,additionals);
subj_stab_all=subj_stab_all(:,2:end);

% just in case additionals are set to true
jc=jc(:,2:end);     
jc(find(isnan(jc)))=0;
jd=jd(:,2:end);     
jd(find(isnan(jd)))=0;
ri=ri(:,2:end);     
ri(find(isnan(ri)))=0;
vi=vi(:,2:end);     
vi(find(isnan(vi)))=0;
nmi=nmi(:,2:end);   
nmi(find(isnan(nmi)))=0;



% plot the results

fig= figure;
set(fig, 'Position', [0 0 1300 2500])

subplot(2,3,1);
plot(x,y,'k*')
title('Artificial Data','FontSize',10)
set(gcf, 'position',pos)
set(gca,'xlim',[0 1],'ylim',[0 1])

subplot(2,2,1);
errorbar(krange,sMean,sStd,'o-');
title(['Stability by Clusters'],'FontSize', 10);
xlabel(['# of clusters'],'FontSize', 10)
ylabel(['Stability Value'],'FontSize', 10)
set(gca,'XTick',krange)
set(gca,'XTickLabel',krange)

subplot(2,2,2);
qualVal(qualVal>1e6)=0;
qualStd(qualStd>1e6)=0;
errorbar(krange,qualVal,qualStd,'o-');
title(['Quality by Clusters'],'FontSize',10)
ylabel(['Quality Value'],'FontSize',10)
xlabel(['# of clusters'],'FontSize',10)
set(gca,'XTick',[krange])
set(gca,'XTickLabel',[krange]);

subplot(2,2,3)
imagesc(sortrows(subj_stab_all))
xlabel(['# of clusters'],'FontSize', 10)
title(['Individual Stability by Clusters'],'FontSize',10)
set(gca,'XTick',1:length(krange))
set(gca,'XTickLabel',[krange]);

subplot(2,2,4)
errorbar(mean(subj_stab_all),std(subj_stab_all),'o-');
title(['Mean Individual Stability by Clusters'],'FontSize',10)
xlabel(['# of clusters'],'FontSize', 10)
ylabel(['Mean Stability Value'],'FontSize', 10)
set(gca,'XTick',1:length(krange))
set(gca,'XTickLabel',[krange]);


subplot(2,3,6)
plot(krange,[jc(1,:);jd(1,:);ri(1,:);vi(1,:);nmi(1,:)],'-o');
xlabel(['# of clusters'],'FontSize', 10)
ylabel(['Similarity Metric Value'],'FontSize', 10)
set(gca,'XTick',[krange])
set(gca,'XTickLabel',[krange]);
legend('JC','JD',...
    'RI','VI',...
    'NMI','Location','SouthEastOutside')
title('Other Similarity Metrics','FontSize',10)

% save fig
ide=strrep(num2str(now),'.','');
print(fig,'-djpeg',[pwd '/output/' ide])
save([pwd '/output/' ide '.mat'],'artdata')

