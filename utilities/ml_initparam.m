function param = ml_initparam(param,paramdefault,paramrequired,parents)
%TZ_INITPARAM Initialize parameters.

%   PARAM2 = TZ_INITPARAM(PARAM,PARAMDEFAULT,PARAMREQUIRED,PARENTS)
%   returns the initialized parameter. PARAM and PARAMDEFAULT are both
%   structures. All fields in PARAM will be kept unchanged in PARAM2. And
%   all fields in PARAMDEFAULT but not in PARAM will be added to PARAM2.
%   If a field is a structure, parameter initialization will go down
%   hierachically. PARAMREQUIRED is an optional struct with empty fields.
%   If any of these fields is not a field in PARAM, an error is thrown.
%   PARENTS is for internal use.

%   
%   Example:
%       param2 =  ml_initparam(struct('t',1,'t2',2),struct('t2',3,'t3',4))
%       returns
%           param2 = 
%               t: 1
%               t2: 2
%               t3: 4
%     
%       param2 =  ml_initparam(struct('t',1,'t2',struct('t3',2)), ...
%               struct('t1',3,'t2',struct('t4',4)))
%       returns
%           param2 = 
%               t: 1
%               t2: struct('t3',2,'t4',4)
%               t1: 3
%
%   See also

% Copyright (C) 2006, 2021  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%   18-Nov-2005 Initial write T. Zhao
%   21-Jul-2021 Added paramrequired and parents arguments T. Buck
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isempty(param)
    param = paramdefault;
    return;
end

if ~exist('paramrequired', 'var')
    paramrequired = struct();
end

if ~exist('parents', 'var')
    parents = '';
end

requiredParameterNames = fieldnames(paramrequired);

for k=1:length(requiredParameterNames)
    name = requiredParameterNames{k};
    full_name = name;
    if ~isempty(parents)
        full_name = [parents, '.', full_name];
    end
    if ~isfield(param,name)
        error('%s is required', full_name);
    end
end

defaultParameterNames = fieldnames(paramdefault);

for k=1:length(defaultParameterNames)
    name = defaultParameterNames{k};
    defaultParameterValue = getfield(paramdefault,name);
    if ~isfield(param,name)
        param = setfield(param,name,defaultParameterValue);
    else
        if isstruct(defaultParameterValue)
            parameterValue = getfield(param,name);
            paramrequired2 = struct();
            if isfield(paramrequired,name)
                paramrequired2 = getfield(paramrequired,name);
            end
            parents2 = name;
            if ~isempty(parents)
                parents2 = [parents, '.', parents2];
            end
            param = setfield(param,name,ml_initparam( ...
                parameterValue,defaultParameterValue, ...
                paramrequired2, parents2));
        end
    end
end