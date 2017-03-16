function [sMean,sStd,subj_stab_all]=eval_clustrs(X, krange, dist, runs)
disp(['Cluster Stability, Benhur method...'])
S=cluster_stability(X,krange,0.8,runs,dist);
S(1,:)=[]; sMean=mean(S,2); sStd=std(S,0,2);


% subject stability
disp(['calculating subject stability...'])
IDX=zeros(size(X,1),runs,length(krange));
for k=krange,
    for m=1:runs,
        IDX(:,m,k) = kmeans(X, k,'Replicates',5,'distance',dist);
        disp([num2str(m) ' - ' num2str(k)])
    end
end


subj_stab_all=zeros(size(X,1),max(krange)-1);
for k=krange,
    cube=zeros(size(X,1),size(X,1),runs);
    for r=1:runs
        PART1=IDX(:,r,k-1); %IDX(:,r);
        C1=zeros(length(PART1),length(PART1));
        for i1=1:length(PART1),
            for j1=1:length(PART1),
                if PART1(i1,1)==PART1(j1,1),
                    C1(i1,j1)=1;
                end
            end
        end
        cube(:,:,r)=C1;
    end
    
    subj_stab=zeros(size(X,1),1);
    for s=1:size(X,1),
        crosscorr=corrcoef(reshape(cube(s,:,:),size(X,1),runs));        
        upper_tri=triu(crosscorr);
        upper_tri_vec=crosscorr(find(upper_tri(:)~=0))';
        subj_stab(s) = mean(upper_tri_vec);
        
    end
    disp(k);
    subj_stab_all(:,k-1) = subj_stab;
end
