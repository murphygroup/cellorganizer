function [ cellparam ] = img2param(imdna, imcell, improt, mask, savedir, options)
%IMG2PARAM Turns an image into a parameterization for Cellorganizer

% Gregory Johnson
%
% Copyright (C) 2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

% icaoberg 2016/04/21 Updated method to include the running time in the
% parameterization.

if ~exist('options', 'var')
    options = [];
end

%inital options struct
options = ml_initparam(options, struct('model', []));

if ~exist(savedir, 'dir')
    mkdir(savedir)
end

%set the default dimensionality to the dimensionality of the data
if ~isfield(options.model, 'dimensionality')
    if size(imdna,3) == 1
        options.model.dimensionality = '2D';
    else
        options.model.dimensionality = '3D';
    end
end

switch options.dimensionality
    case '2D'
        %do 2D preprocessing
        cellparam = img2param_2D(imdna, imcell, improt, mask, savedir, options);
    case '3D'
        %do 3D preprocessing
        tic;
        cellparam = img2param_3D(imdna, imcell, improt, mask, savedir, options);
        
        %this is done so we can measure how long did it take to measure the
        %parameterization of the cell
        if ~isempty(cellparam)
            cellparam.options.running_time.parameterization = toc;
        end
    otherwise
        disp(['Unsupported dimensionality ' options.model.dimensionality '. Returning empty parameterization.'])
end
end