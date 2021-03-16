function ts = tz_expandfunname(funname)
%TZ_EXPANDFUNNAME Expand function name
%   TS = TZ_EXPANDFUNNAME(S) returns expansion of FUNNAME.
%   It is useful for taking a quick look at how the function or a variable
%   is named :).
%   TS is cell array.

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

funnameTokens = tz_strtok(funname,'_');
ts = {};

for i = 1:length(funnameTokens)
    fullWord = tz_expandcombabbrv_fast(funnameTokens{i});
    ts = {ts{1:end},fullWord{1:end}};    
end

disp( tz_cell2str(ts,' ') );
