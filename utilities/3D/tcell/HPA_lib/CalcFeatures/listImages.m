function [readlist] = listImages(rootdir, writedir, naming_convention)

% rootdir is the root directory storing all images.

% Copyright (C) 2010  Murphy Lab
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

% 2011-11-12 tebuck: copied from calcRegionFeat.m just to get the
% list of images considered (just channel zero)


fsetnames = {...
    'overlap', ...
    'overlapx2', ...
    'nonObjFluor', ...
    'nonObjFluorx2', ...
    'obj', ...
    'objx2', ...
    'mutualInfo', ...
    'mutualInfox2', ...
    'texture', ...
    'texturex2', ...
    'texturex4', ...
    'tas', ...
    'tasx2', ...
    'objRegion',...
    'objRegionx2',...
};


datasettype = 'Region';
% fsetnames = {'nuclearRegion'};
% $$$ if ~exist(['./data/features/' lower(datasettype) '/'],'dir')
% $$$    mkdir(['./data/features/' lower(datasettype) '/']);
% $$$ end
if ~exist([writedir 'features/' lower(datasettype) '/'],'dir')
   mkdir([writedir 'features/' lower(datasettype) '/']);
end


if ~exist('rootdir','var')
    rootdir = '/images/HPA/images/IFconfocal/';
end
if ~exist('naming_convention', 'var') || strcmpi(naming_convention, 'IFconfocal')
  naming_convention = struct(); 
  naming_convention.protein_channel = 'green'; 
  naming_convention.nuclear_channel = 'blue'; 
  naming_convention.tubulin_channel = 'red'; 
  naming_convention.er_channel = 'yellow'; 
elseif exist('naming_convention', 'var') && strcmpi(naming_convention, 'IFconfocal_CellCycle')
  naming_convention = struct(); 
  naming_convention.protein_channel = '_ch00'; 
  naming_convention.nuclear_channel = '_ch01'; 
  naming_convention.tubulin_channel = '_ch03'; 
  naming_convention.er_channel = '_ch02'; 
end

if ~isfield(naming_convention, 'segmentation_suffix')
  naming_convention.segmentation_suffix = naming_convention.protein_channel; 
end

filetype = 'tif';
%greparg = '| grep green';
%greparg = ['| grep ', naming_convention.protein_channel];
greparg = ['| grep ', naming_convention.segmentation_suffix];

ind = find(rootdir=='/');
rootdir_ = rootdir;
rootdir_(ind) = '_';

uout = unixfind( rootdir, filetype, greparg);
readlist = listmatlabformat( uout);

