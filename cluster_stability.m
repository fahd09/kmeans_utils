function [S]=cluster_stability(CMAT,krange,prop,simreps,dist)
%cluster_stability     Cluster Stability
%
%   S = cluster_stability(CMAT,krange,prop,simreps);
%   This function calculates the cluster stability across range of K's.
%
%   Inputs: CMAT,  MxN matrix where M is # of observations (rows) and N
%   is # of variables.
%           krange,  1xP vector of the values which stability is estimated
%           prop,  proportion of dataset to be used in each repetition
%           simreps,  The nunmber of times stability estimated at each k
%
%   Output: S,  KxR matrix where K represents cluster k and R is
%   the number of repetitions.
%
%   Example:
%            krange=[2:10]; prop=0.8;simreps=100;
%            S=cluster_stability(DATA,krange,prop,simreps);
%            errorbar(krange,mean(S,2),std(S,0,2),'-o');
%
% Reference:
% Salvador S and Chan P. Determining the number of clusters/segments in hierarchical clustering/
% segmentation algorithms. November 2004.
% Ben-Hur A, Elisseeff A, and Guyon I. A stability based method for discovering structure in
% clustered data. Pacific Symposium on Biocomputing, 7:6–17, 2002.


num_samps=round(prop*size(CMAT,1));
S=zeros(length(krange),simreps);

for k=krange,
    
    for m=1:simreps,
        
        % generating sub-sample indices
        
        r1=randperm(size(CMAT,1))';
        ss1=r1(1:num_samps,1);
        r2=randperm(size(CMAT,1))';
        ss2=r2(1:num_samps,1);
        
        [idx1]=kmeans(CMAT(ss1,:),k,'Start','sample','Replicates',5,'distance',dist);
        [idx2]=kmeans(CMAT(ss2,:),k,'Start','sample','Replicates',5,'distance',dist);
        
        
        % find set intersection
        
        [~,iss1,iss2]=intersect(ss1,ss2);
        
        D1=idx1(iss1,1);
        D2=idx2(iss2,1);
        
        %% computing similarity between matrices
                
        S(k,m)=partsim(D1,D2);
        disp([num2str(m) ' - ' num2str(k)])
        
        
    end
    
    % displaying
    
    % display(k);
    
end

end