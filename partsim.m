function [PS]=partsim(PART1,PART2)
%partsim     Partition Similarity
%
%   PS = partsim(PART1,PART2);
%   This function finds the similarity between two cluster assignments.
%
%   Inputs: PART1 & PART2,  Two vectors of community assignments.
%   Output: PS,  the correlation between the two assignments.
%

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
