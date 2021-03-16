function [ h, x, y ] = plotmodelsparam( models, x, modelvar, plotfunc )
%PLOTMODELSPARAM plots the independent variable, x, against the models struct component, modelVar.
%The model parameter value gets concatonated along the first non-singleton
%dimension
%
%Inputs:
% models = a cell array of models to compare
% x = the independent variable
% modelvar = the field of the model to compare
%plotfunc = (optional) function handle for plotting the parameters
%
%Outputs:
% h = plot handle
% x = array of x values
% y = array of y values

%Author: Greg Johnson 16/5/2013
%Edited: D. Sullivan 9/27/13 - documentation
%
% Copyright (C) 2013 Murphy Lab
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

if ~exist('x', 'var')
    x = ones(1, length(models));
end

if iscell(models)
    models = [models{:}];
end

tokens = regexp(modelvar, '\.', 'split');

%Get the data from each model
for i = 1:length(models)
    model = models(i);
    for j = 1:length(tokens)
        model = model.(tokens{j});
    end
   
    y{i} = model;
end

%Concationate data by the first singleton dimension
ysize = size(y{1});

inds = find(ysize ==1);

if isempty(inds)
    inds = length(ysize) +1;
end

y = cat(inds(1), y{:});
        
if exist('plotfunc', 'var') & ~isempty(plotfunc)
    h = plotfunc(x, y);
else
    h = scatter(x,y);
end

end
