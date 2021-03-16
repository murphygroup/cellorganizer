function ids = tz_genmixsel(n,s,m,rns,pmin)
%TZ_GENMIXSEL Randomly generate mixture indices.
%   IDS = TZ_GENMIXSEL(N,S,M) returns the cell array of random generated
%   indices of cells. N is the total numbers of cells in the patterns. S 
%   is the selected numbers of cells in the patterns. N and S must have
%   the same length except the case when S is empty, in which 1 to 5 cells
%   will be selected for each pattern. If S(I) is negative, then 1 to 5 
%   cells will be selected for pattern I. IDS is a cell array with each
%   element for each pattern. Each element is a vector in which each 
%   number is the index of a cell except 0, which means nothing.
%   
%   IDS = TZ_GENMIXSEL(N,S,M,RNS) lets user specify the possible number of
%   cells in certain patterns. RNS is only effective for negative elements
%   in S or for all patterns if S is empty. TZ_GENMIXSEL(N,S,M) is the same
%   as TZ_GENMIXSEL(N,S,M,1:5).
%   
%   IDS = TZ_GENMIXSEL(N,S,M,RNS,PMIN) allow setting the minimal number of
%   patterns, which will not be less than PMIN.

%   17-May-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University
  
if nargin < 3
    error('At least 3 arguments are required')
end

if ~exist('pmin','var')
    pmin=1;
end

if pmin>length(n)
    error('minimal number of patterns must not be greater total number of patterns');
end

if ~exist('rns','var')
    rns=1:5;
end

if any(rns<0)
    error('rns must not be less than 0');
end

if any(max(rns)>n)
    error('rns must not be greater than total number of cells');
end

if isempty(s)
    s=zeros(size(n))-1;
end

if sum(s~=0)<2
    error('At least two samples should be selected');
end

rs=s;

if all(s==1)
    tmixsize=prod(n);
    mixperm=randperm(tmixsize);
    ss=[];
    for i=1:length(n)
        ss=[ss ' ids{' num2str(i) '},']; 
    end
    eval(['[' ss(1:end-1) ']=ind2sub(n,mixperm(1:m)'');']);

    return
end

%go through each pattern
for i=1:length(n)
    if rs(i)<0  %generate random number of cells
        for j=1:m
            tmps=randperm(length(rns));
            s(i)=rns(tmps(1));
            if s(i)==0
                ids{i}(j,1)=0;
            else
                mixperm=randperm(n(i));
                for k=1:s(i)
                    ids{i}(j,k)=mixperm(k);
                end
            end
        end
    else
        if s(i)==1
            mixperm=randperm(n(i));
            nm=length(mixperm);
            ids{i}=mixperm(mod(1:m,nm)+1)';
        else
            if s(i)==0
                ids{i}=zeros(m,1);
            else
                for j=1:m
                    mixperm=randperm(n(i));
                    for k=1:s(i)
                        ids{i}(j,k)=mixperm(k);
                    end
                end
            end
        end
    end
end

for i=1:length(ids)
    if size(ids{i},2)<max(rns)
        ids{i}=[ids{i},zeros(size(ids{i},1),max(rns)-size(ids{i},2))];
    end
end

combids=[];
for i=1:length(ids)
    combids=[combids,sum(ids{i},2)>0];
end

%iteration. replace ids unsatisfying with the conditions
invalids=find(sum(combids,2)<pmin);
if ~isempty(invalids)
    tmpids = tz_genmixsel(n,rs,length(invalids),rns,pmin); 
    for i=1:length(ids)
        if sum(sum(tmpids{i},1)>0)>0   
%             sum(tmpids{i},1)>0
            ids{i}(invalids,1:size(tmpids{i},2))=tmpids{i};
            if size(ids{i},2)>size(tmpids{i},2)
                ids{i}(invalids,size(tmpids{i},2)+1:end)=0;
            end
        end
    end
end
