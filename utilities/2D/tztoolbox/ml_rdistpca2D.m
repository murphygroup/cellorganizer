function [f,avgratio] = ml_rdistpca2D(combcellcodes,param)
%ML_RDISTPCA Train the the PCA model for distance ratios.
%   F = ML_RDISTPCA(COMBCELLCODES) returns the trained PCA model for the
%   input cell array of cell codes COMBCELLCODES. F is a structure with the
%   following fields:
%       'startangle' - see below.
%       'anglestep' - see below.
%       'minratio' - minimal ratios. It will be ignored if it is empty. If it
%          is 0, the minimal ratio will be found automatically.
%       'stat' - a pdf.
%   
%   F = ML_RDISTPCA(COMBCELLCODES,PARAM) specifies parameters for training.
%   PARAM is a structure that has the following fields:
%       'startangle' - a string that determines how to align the ratio
%           vector. 'cell' means the major angle of the cell and 'nuc'
%           means the major angle of the nucleus.
%       'anglestep' - step of the angles. It must be an integer.
%       'minratio' - if it less than 1, the ratio will be obtained from the
%           data.
%       'ml_estpdf' - parameters for the function ML_ESTPDF
%
%   [F,AVGRATIO] = ML_RDISTPCA(...) also returns the average shape ratio.
%
%   See also

% Copyright (C) 2007-2013 Murphy Lab
% Center for Bioimage Informatics
% Lane Center for Computational Biology
% Carnegie Mellon University
%
% 10-Jan-2006 Initial write T. Zhao
% 07-Mar-2013 I. Cao-Berg Updated method so that if the distance vector in the
%             cell codes structure is not of size 360, it ignores that image for
%             training
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

if nargin < 1
    error('Exactly 1 argument is required')
end

if ~exist('param','var')
    param = struct([]);
end

%icaoberg 5/10/2013
if ~isfield(param, 'anglestep' )
%     if length(combcellcodes)>360
        param.anglestep=1;
%     else
%         param.anglestep = round(359/length(combcellcodes)+0.5);
%     end
end

param = ml_initparam(param, ...
    struct('startangle','nuc','minratio',[],'ml_estpdf',struct([])));
c = 1;

%grj 14/5/2013
%preallocate normrdist matrix
try
    cc = [combcellcodes{:}];
%     cc = [cc.cellcode];
    normrdist = zeros(length(cc), median(cellfun(@length, {cc.nucelldist})));
catch
    normrdist = [];
end

for i=1:length(combcellcodes)
    cellcode = combcellcodes{i};
    
    %icaoberg 07/03/2013
    %ignore cell code if distance vectors are not the right size, i.e. 1x360
    %grj      14/05/2013
    %added AND to make sure everything fits in normrdist
    %added AND to make sure rdist doesnt go to inf and break ml_estpdf2D
    if size( cellcode.nucelldist, 2 ) == size( cellcode.nucdist, 2 ) ...
            && length(cellcode.nucelldist) == size(normrdist(c,:),2) ...
            && ~(any(cellcode.nucdist == 0))
        rdist = cellcode.nucelldist./cellcode.nucdist;
        
        switch param.startangle
            case 'cell'
                startAngle = cellcode.cellmangle;
            case 'nuc'
                startAngle = cellcode.nucmangle;
        end
        %rdist = rdist((0:param.anglestep:359)+1);
        
        if ~isnan(startAngle)
            normrdist(c,:) = ml_shiftdist(rdist,startAngle);
        else
            continue;
        end
        
        c = c+1;
    end
end

%grj 14/05/2013
%remove the empty cells of normrdist
normrdist(c:end,:) = [];

f.startangle = param.startangle;
f.anglestep = param.anglestep;
if ~isempty(param.minratio)
    if param.minratio<1
        f.minratio = min(normrdist(:));
    else
        f.minratio = param.minratio;
    end
else
    model.minratio = [];
end

avgratio = mean(normrdist,1);


%normrdist = ml_addrow(normrdist,-avgratio);
%f.avgratio = avgratio;

% if isfield(param.ml_estpdf.transform.param,'ncomp')
%     param.ml_estpdf.mu = zeros(1,param.ml_estpdf.transform.param.ncomp);
% else
%     param.ml_estpdf.mu = zeros(1,size(normrdist,2));
% end

f.stat = ml_estpdf2D(normrdist,param.ml_estpdf);

