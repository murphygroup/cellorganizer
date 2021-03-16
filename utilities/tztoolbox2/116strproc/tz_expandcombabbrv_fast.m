function ts = tz_expandcombabbrv(funname)
%TZ_EXPANDCOMBABBRV_FAST Expand combination of abbreviations by a fast alg.
%   TS = TZ_EXPANDCOMBABBRV_FAST(S) returns full words of S, which is the
%   combination of a set abbreviations. It is useful for taking a quick look 
%   at how the function or a variable is named.

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

gm = tz_buildabbrvgm(funname);
nodes = ml_shortpath(gm);

fullword2 = [];
k = 1;
flag = 0;

for i=1:length(nodes)-1
    subs = funname(nodes(i):nodes(i+1)-1);
    [fullword,stat] = tz_expandabbrv(subs);
    
    if stat==1
        if ~isempty(fullword2)
            ts{k} = fullword2;
            fullword2 = [];
            k = k+1;
        end
        flag = 1;
        ts{k} = fullword;
        k = k+1;
    else
        fullword2 = [fullword2 subs];
        flag = 0;
    end
end

if flag==0
    ts{k} = fullword2;
end



