function img_out = param2img( param, options )
%PARAM2IMG Returns an instance from the Parametrization  of a cell. It
%reads the parameters in the variable 'parametrization' and generates an
%instance

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: October 12, 2015
%
% Copyright (C) 2015 Murphy Lab
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

%% step0: check input arguments
%check number of input arguments

if ~exist('options', 'var')
    options = [];
end

% xruan 01/05/2016 check if the param is the filename or the struct
if ischar(param) && exist(param, 'file')
    param = load(param);
end

if ~isstruct(param)
    error(['the percell parameter must be a struct array or its filename']);
end

options = ml_initparam(options, struct('components', {fieldnames(param)}));

for i = 1:length(options.components)
    %the components parameters come first, then the parameters as follows
    %in numerical order in param.dependencies
    
    img{i} = synth_img(param.(options.components{i}));
end

% current only implement indexed image

switch(options.format)
    case 'indexed'
        
        for i = 1:length(options.components)
            img_i_out = img{i};
            if isfield(param.options.model, 'downsampling')
                img_i_out = ml_downsize(img_i_out, param.options.model.downsampling, 'linear');
                img_vals = unique(img{i}(:));
                
                if numel(img_vals) == 2
                    img_i_out = img_i_out >= (img_vals(img_vals > 0) / 2);
                end                    
            end
            if i == 1 
                img_out = img_i_out;
            else                
                img_out = img_out + img_i_out;
            end
        end

    otherwise
        error('only implement indexed case, rgb, greyscale etc will be implement later!')
end

end%param2img

function component_img = synth_img(component_param)

    switch(component_param.type)
        case 'medial axis'
        case 'cylindrical surface'
        case 'ratio'
        case 'gaussian'
        case 'network'
        
        % xruan add case for diffeomorphic, currently only for
        % training. latter and synthesis. 
        case 'diffeomorphic'
            if isfield(component_param, 'seg')
                component_img = component_param.seg;
            else
                error('synthesis not implement yet!')
            end
        otherwise 
            error(['Unrecognized component type ' component_param.type '.'])
    end

end