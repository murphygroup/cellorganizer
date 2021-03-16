function [ model ] = param2model(cell_params, options)
% PARAM2MODEL Turns a bunch of per-cell-parameterizations into a model of a distribution
% over these parameters

% Copyright (C) 2012-2019 Murphy Lab
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

if ~exist('options', 'var')
    options = [];
end

%inital options struct
options = ml_initparam(options, struct());

if ~isfield(options,'dimensionality')
    error('No dimensionality specified')
end

switch options.dimensionality
    case '2D'
        %do 2D preprocessing
        model = param2model_2D(cell_params, options);
        model.dimensionality = options.dimensionality;
    case '3D'
        %do 3D preprocessing
        model = param2model_3D(cell_params, options);
        model.dimensionality = options.dimensionality;
    otherwise
        disp(['Unsupported dimensionality ' options.dimensionality '. Returning empty model.'])
        model = [];
end

if isfield(options,'dataset')
    model.dataset = options.dataset;
end

if isempty(model)
    warning('Model is empty.');
    return
end
end%param2model.m
