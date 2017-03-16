function [jaccardcoeff jaccarddistance randindex variationofinformation normalizedmutualinformation] = similarity_indices(clrfile1,clrfile2,varargin)

% jdp 10/10/10
% 
% This script takes two identically sized sets of module assignments, and
% computes common similarity indices on them, including:
% 
% jaccard coefficient (JC)
% jaccard distance (JD)
% rand index (RI)
% variation of information (VI)
% normalized mutual information (NMI)
% 
% These are common measures; google will familiarize the uninitiate.
% 
% clrs can be either files or arrays, and passing in a 3rd argument
% suppresses written output.
% 
% [jc jd ri vi nmi] = similarity_indices(clrfile1,clrfile2)
% [jc jd ri vi nmi] = similarity_indices('clrfile1.txt','clrfile2.txt')
% [jc jd ri vi nmi] = similarity_indices(clrmat1,clrmat2,0)

clf;
if isempty(varargin)
    writeoutput=0;
else
    writeoutput=1;
end

% load clrfiles if necessary
clrs1=clrfile1;
clrs2=clrfile2;
if ~isnumeric(clrfile1)
    clrs1=load(clrfile1);
end
if ~isnumeric(clrfile2)
    clrs2=load(clrfile2);
end

% ensure dimensions match
d1 = size(clrs1);
d2 = size(clrs2);
if ~isequal(d1,d2)
    error('Clrfiles not same size');
end

% for each column calculate rand and jaccard indices, and infotheory stuff
for i=1:d1(2)
    [jaccardcoeff(i,1) jaccarddistance(i,1) randindex(i,1) ] = jaccard_rand_indices(clrs1(:,i),clrs2(:,i));
    [variationofinformation(i,1) normalizedmutualinformation(i,1)] = ITvi_nmi(clrs1(:,i),clrs2(:,i));
end

% write output if output isn't suppressed with varargin
if writeoutput
    
    % create filestem for the output
    [pathstr,name1,ext,versn] = fileparts(clrfile1);
    [pathstr,name2,ext,versn] = fileparts(clrfile2);
    riname = [name1 '_' name2 '_similarity_indices'];
    
    % save the indices
    fid=fopen([riname '.txt'],'w'); fprintf(fid,'jaccardcoeff\tjaccarddistance\trandindex\tvariationofinformation\tnormalizedmutualinformation\n'); fclose(fid);
    dlmwrite([riname '.txt'],[jaccardcoeff jaccarddistance randindex variationofinformation normalizedmutualinformation],'delimiter','\t','-append');
    
    % plot and save the indices
    range=[1:d1(2)];
    subplot(1,2,1);
    plot(range,jaccardcoeff(:,1),'b-',range,jaccarddistance(:,1),'r-',range,randindex(:,1),'g-');
    xlabel('Thr/box'); ylabel('JC-blue JD-red RI-green');
    subplot(1,2,2);
    plot(range,variationofinformation(:,1),'k-',range,normalizedmutualinformation(:,1),'c-');
    xlabel('Thr/box'); ylabel('VI-black NMI-cyan');
    saveas(gcf,[riname '.tiff'],'tiff');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [jc jd ri] = jaccard_rand_indices(col1,col2)

[tallymat1] = bootstrapper_tally(col1);
[tallymat2] = bootstrapper_tally(col2);

tallymat1=triu(tallymat1,1);
tallymat2=triu(tallymat2,1);

[a b]=size(tallymat1);
tot=((a*a)-a)/2;

M10=nnz(find((tallymat1-tallymat2)>0));
M01=nnz(find((tallymat1-tallymat2)<0));
M11=nnz(find((tallymat1+tallymat2)>1));
M00=tot-M11-M10-M01;

jc=M11/(M01+M10+M11);
jd=(M10+M01)/(M01+M10+M11);
ri=(M11+M00)/(M00+M01+M10+M11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vi nmi] = ITvi_nmi(col1,col2)

nodes1=size(col1,1);
nodes2=size(col2,1);

% find joint probability
[p12 p1 p2] = jointprobability(col1,col2);

% calculate entropies
h1=ITentropy(p1);
h2=ITentropy(p2);

% get the mutual information
[mi] = ITmutualinfo(p12,p1,p2);

vi=h1+h2-(2*mi);

vi=vi/(log2(nodes1));

nmi = (2*mi)/(h1+h2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mi] = ITmutualinfo(p12,p1,p2)

[a1 b1]=size(p1);
[a2 b2]=size(p2);
[a12 b12]=size(p12);

if ~(isequal(a1,a12) && isequal(a2,b12))
    error('Trouble, probability matrices don''t seem the right sizes');
end

mi=0;
for i=1:a12
    for j=1:b12
        if p12(i,j)~=0 % you're multiplying Inf * 0, which is 0
            mi=mi+(p12(i,j)*log2(p12(i,j)/((p1(i,1)*p2(j,1)))));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [entropy] = ITentropy(p)

% how many states?
[a b]=size(p);

entropy=0;
for i=1:a
    entropy=entropy+(p(i)*log2(p(i)));
end

entropy=entropy*(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p] = column2probabilitydistribution(column)

nodes=size(column,1);
mods=unique(column);
nummods=size(mods,1);
p=zeros(nummods,1);
for m=1:nummods
    p(m,1)=nnz(find(column==mods(m)));
end
p=p/nodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p12 p1 p2] = jointprobability(col1,col2)

nodes1=size(col1,1);
nodes2=size(col2,1);

mods1=unique(col1);
nummods1=size(mods1,1);
mods2=unique(col2);
nummods2=size(mods2,1);

p1=column2probabilitydistribution(col1);
p2=column2probabilitydistribution(col2);

for i=1:nummods1
    for j=1:nummods2
        n12(i,j)=nnz(double(ismember(find(col1==mods1(i)),find(col2==mods2(j)))));
    end
end

% define joint probability
p12=n12/nodes1;


