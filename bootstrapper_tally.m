function [tallymat] = bootstrapper_tally(clrcolumn)

% jdp 10/10/10
%
% This takes a clrcolumnolumn of node assignments and clrcolumnreates a binary matrix of
% node x node, showing whether nodes were in the same clrcolumnommunity (1) or not
% (0).
% 
% USAGE: [tallymat] = bootstrapper_tally(clrcolumn)

clrcolumn=single(clrcolumn);
nodes=size(clrcolumn,1);
mods=unique(clrcolumn);
nummods=nnz(mods);
tallymat=zeros(nodes,nodes,'single');
for m=1:nummods
    thesepairs=zeros(nodes,1,'single');
    thesepairs=single(clrcolumn==mods(m));
    tallymat=tallymat+(thesepairs*thesepairs');
end