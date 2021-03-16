function tz_save(filepath,vars,varnames,script,comments,options)
%TZ_SAVE Save variables with related information
%   TZ_SAVE(FILEPATH,VARS,VARNAMES,SCRIPT,COMMENTS) save the variable in 
%   the cell  array VARS in the file FILEPATH. The names of the variables 
%   are in the [string array] VARNAMES. SCRIPT is the script which created
%   the variables. COMMENTS is a string of comments for the data.
%   
%   TZ_SAVE(FILEPATH,VARS,VARNAMES,SCRIPT,COMMENTS,OPTIONS) also specifies 
%   options allowed in the the function SAVE. OPTIONS is a [string array].
%   
%   See also

%   24-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 5
    error('5 or 6 arguments are required')
end

if isempty(which(script))
    warning(['The script ' script ' not found. Please check if the name ' ...
                        'is correct']);
end

for i=1:length(varnames)
    eval([varnames{i} '= vars{i};']);
end

machine = tz_machineinfo('name','computer','domain');
savetime = clock;
matlab_version = version;

if ~exist('options','var')
    options = {};
end

save(filepath,varnames{:}, ...
    'script','comments','machine','savetime','matlab_version',...
    options{:});
