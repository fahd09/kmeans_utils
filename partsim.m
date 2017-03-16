function [PS]=partsim(PART1,PART2)
%partsim     Partition Similarity
%
%   PS = partsim(PART1,PART2);
%   This function finds the similarity between two cluster assignments.
%
%   Inputs: PART1 & PART2,  Two vectors of community assignments.
%   Output: PS,  the correlation between the two assignments.
%
% Reference:
% Ben-Hur A, Elisseeff A, and Guyon I. A stability based method for discovering structure in
% clustered data. Pacific Symposium on Biocomputing, 7:6–17, 2002.


% computing similarity between matrices
% there are better ways of doing this.
    
% cluster 1
    
%C1=zeros(length(PART1),length(PART1));
% for i1=1:length(PART1),
%         
%      for j1=1:length(PART1),
%             
%          if PART1(i1,1)==PART1(j1,1),
%             
%              C1(i1,j1)=1;
%              
%          end
%             
%     end
%         
% end
    
% cluster 2
    
%C2=zeros(length(PART2),length(PART2));
 
% for i2=1:length(PART2),
%         
%      for j2=1:length(PART2),
%             
%          if PART2(i2,1)==PART2(j2,1),
%             
%               C2(i2,j2)=1;
%              
%          end     
%             
%      end
%         
% end


C1=bootstrapper_tally(PART1);
C2=bootstrapper_tally(PART2);

% dot-product
    
C12=C1.*C2;
    
% similarity terms
    
L12=sum(sum(C12));
L11=sum(sum(C1));
L22=sum(sum(C2));

% partition similarity
    
PS=L12/sqrt(L11*L22);

end