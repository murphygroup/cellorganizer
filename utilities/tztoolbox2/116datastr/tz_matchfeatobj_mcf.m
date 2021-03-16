function succ=tz_matchfeatobj_mcf(all_features,objs)

%TZ_MATCHFEATOBJ_MCF Check if cell features and object are matched.
%   SUCC = TZ_MATCHFEATOBJ_MCF(ALL_FEATURES,OBJS) returns 1 if 
%   ALL_FEATURES has can be matched OBJS. Otherwise 0 will be
%   return. ALL_FEATURES are MCF and the first column in each matrix
%   should contain object numbers. OBJS is the three-level cell array
%   of objects.

%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZS
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if nargin < 2
    error('Exactly 2 arguments are required')
end

nclass=length(all_features);

if nclass~=length(objs)
    succ=0;
    return
end

for i=1:nclass
    ncell=size(all_features{i},1);
    if ncell~=length(objs{i});
        succ=0;
        return
    end
    
    for j=1:ncell
        if all_features{i}(j,1)~=length(objs{i}{j})
            succ=0;
            warning(['class ' num2str(i) 'cell ' num2str(j) ' unmatched']);
            [all_features{i}(j,1) length(objs{i}{j})]
            %return;
        end
    end
end

succ=1;