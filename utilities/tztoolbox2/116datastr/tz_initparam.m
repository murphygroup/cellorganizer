function param = tz_initparam(param,paramdefault)
%TZ_INITPARAM Initialize parameters. (Obsolete)
%
%See also ML_INITPARAM

%   PARAM = TZ_INITPARAM(PARAM,PARAMDEFAULT) returns the initialized 
%   parameter. PARAM and PARAMDEFAULT are both structures.
%   
%   See also

%   18-Nov-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_initparam','ml_initparam'));

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isempty(param)
    param = paramdefault;
    return;
end

defaultParameterNames = fieldnames(paramdefault);

for k=1:length(defaultParameterNames)
    if ~isfield(param,defaultParameterNames{k})
        defaultParameterValue = ...
            getfield(paramdefault,defaultParameterNames{k});
        param = setfield(param,defaultParameterNames{k}, ...
            defaultParameterValue);
    end
end