function ids = tz_genmixsel(n,s,m,smin,smax,pmin)
%TZ_GENMIXSEL Obsolete.

%function s = tz_genmixsel(n,s,m,smin,smax,pmin)
%OVERVIEW
%   generate mixture indices
%PARAMETERS
%   n - totals number in the patterns
%   s - selected numbers in the patterns
%   m - number of generated samples
%RETURN
%   ids - cell array of indices
%DESCRIPTION
%   
%HISTORY
%   17-May-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_genmixsel_bk','tz_genmixsel'));

if ~exist('pmin')
    pmin=1;
end

if pmin>length(n)
    error('minimal number of patterns must not be greater that total number of patterns');
end

if ~exist('smax','var')
    smin=1;
end

if ~exist('smax','var')
    smax=5;
end

if smin<0
    error('smin must not be less than 0');
end

if any(smax>n)
    error('smax must not be greater than total number of cells');
end

if isempty(s)
    s=zeros(size(n))-smax;
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
            s(i)=tz_unifrnd_discr(1,-rs(i)-smin+1)+smin-1;
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

combids=[];
for i=1:length(ids)
    combids=[combids,sum(ids{i},2)>0];
end

%iteration. replace ids unsatisfying with the conditions
invalids=find(sum(combids,2)<pmin);
if ~isempty(invalids)
    tmpids = tz_genmixsel(n,rs,length(invalids),smin,smax,pmin); 
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