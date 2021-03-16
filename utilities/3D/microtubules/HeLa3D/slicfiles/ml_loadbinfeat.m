function [feat,featids,slfnames,featnames,featdesps,sampleids] = ...
    ml_loadbinfeat(featfile,machineformat)
%ML_LOADBINFEAT  Load binary feature file
%   FEAT = ML_LOADBINFEAT(FILEPATH) load a file from FILEPATH to get a
%   feature matrix FEAT. The file should start with two integers specifying
%   the size of the feature matrix. Each feature should saved as a float
%   in the file.
%   
%   FEAT = ML_LOADBINFEAT(FILEPATH,MACHINEFORMAT) let the user to specify
%   the machine format. See FREAD for more details.
%
%   [FEAT,FEATIDS,SLFNAMES,FEATNAMES,FEATDESPS,SAMPLEIDS] =
%   ML_LOADBINFEAT(:) also returns features IDs, SLF names, feature names,
%   the description of the features and sample IDs.
%
%   See also

% Copyright (C) 2006  Murphy Lab
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

%   16-Oct-2005 Initial write T. Zhao

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('machineformat','var')
    machineformat = 'ieee-be';
end

fid = fopen(featfile,'r');

version = fread(fid,1,'int32',machineformat);
nsample = fread(fid,1,'int32',machineformat);
nfeat = fread(fid,1,'int32',machineformat);
feat=fread(fid,nsample*nfeat,'float32',machineformat);
feat = reshape(feat,nfeat,nsample)';

featids = [];
slfnames = {};
featnames = {};
featdesps = {};
sampleids = [];

switch version
    case {1,2,3,4,5}
        %read feature ids
        if ismember(version,[2 3 4 5])
            nfeatid = fread(fid,1,'int32',machineformat);
            if(nfeatid>0)
                featids = fread(fid,nfeatid,'int32',machineformat);
            end
        end
        
        %read slf names
        breaker = '@';
        infoLength = fread(fid,1,'int32',machineformat);
        if infoLength>0
            slfString = fread(fid,infoLength,'uint8=>char',machineformat);
            slfnames = ml_strtok(slfString',breaker);
        end
        
        %read feature names
        infoLength = fread(fid,1,'int32',machineformat);
        if infoLength>0
            featString = fread(fid,infoLength,'uint8=>char', ...
                               machineformat);
            featnames = ml_strtok(featString',breaker);
        end

        %read feature descriptions
        infoLength = fread(fid,1,'int32',machineformat);
        if infoLength>0
            despString = fread(fid,infoLength,'uint8=>char', ...
                               machineformat);
            featdesps = ml_strtok(despString',breaker);
        end
        
        %read sample Ids
        if ismember(version,[3 4 5])
            nsampleid = fread(fid,1,'int32',machineformat);
            if(nsampleid>0)
                sampleids = fread(fid,nsampleid,'int32',machineformat);
            end
        end
    otherwise
        error('undefined version');
end

fclose(fid);

function ts = ml_strtok(s,d)
%ML_STRTOK Find tokens in string.
%   TS = ML_STRTOK(S) returns all tokens in the string S delimited by
%   "white space". TS is a cell array of these tokens.
%   
%   TS = ML_STRTOK(S,D) returns all tokens delimited by one of the 
%   characters in D. TS is a cell array of these tokens.

if nargin<1
    error('1 or 2 arguments are required');
end

if nargin<2
    d=' ';
end

if isempty(d)
    ts{1}=s;
    return
end

ts={};

[token,remain] = strtok(s,d);

while ~isempty(token)    
    ts{end+1}=token;
    [token,remain] = strtok(remain,d);
end
