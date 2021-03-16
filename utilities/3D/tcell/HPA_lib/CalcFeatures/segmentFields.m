function segmentFields(readdir, writedir, naming_convention)

% function segmentFields
% readdir is the directory storing all images;
% writedir is the directory for writing segmentation masks.

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

% 10 Jan 10 - Jieyue Li
% 2011-11-09 tebuck: adding option to use IFconfocal_CellCycle-like
% image names.
% $$$ % 2011-11-15 tebuck: if naming_convention.nuclear_channel is
% $$$ % 'ignore_channel', it is now set to a blank image the size of the
% $$$ % protein image. Added IFconfocal_CellCycle_nuclear-focus
% $$$ % naming_convention.
% 2011-12-07 tebuck: adding option
% naming_convention.segmentation_suffix to allow reordering
% channels without erroneously affecting segmentation.


if ~exist('readdir','var')
    readdir = '/images/HPA/images/IFconfocal/';
end
if ~exist('writedir','var')
    writedir = './data/masks/';
end
if ~exist('naming_convention', 'var') || strcmpi(naming_convention, 'IFconfocal')
  naming_convention = struct(); 
  naming_convention.protein_channel = 'green'; 
  naming_convention.nuclear_channel = 'blue'; 
  naming_convention.tubulin_channel = 'red'; 
  naming_convention.er_channel = 'yellow'; 
  naming_convention.segmentation_suffix = naming_convention.protein_channel; 
elseif exist('naming_convention', 'var') && strcmpi(naming_convention, 'IFconfocal_CellCycle')
  naming_convention = struct(); 
  naming_convention.protein_channel = '_ch00'; 
  naming_convention.nuclear_channel = '_ch01'; 
  naming_convention.tubulin_channel = '_ch03'; 
  naming_convention.er_channel = '_ch02'; 
  naming_convention.segmentation_suffix = naming_convention.protein_channel; 
end


IMAGEPIXELSIZE = 0.05; % um/px
MINNUCLEUSDIAMETER = 6; %um
MAXNUCLEUSDIAMETER = 14; %um

filetype = 'tif';

%greparg = '| grep green';
greparg = ['| grep ', naming_convention.protein_channel];

ind = find(readdir=='/');
readdir_ = readdir;
readdir_(ind) = '_';

uout = unixfind( readdir, filetype, greparg);
readlist = listmatlabformat( uout);

uout_nuc = findreplacestring(uout, naming_convention.protein_channel, naming_convention.nuclear_channel);
uout_tub = findreplacestring(uout, naming_convention.protein_channel, naming_convention.tubulin_channel);
uout_er = findreplacestring(uout, naming_convention.protein_channel, naming_convention.er_channel);
readlist_nuc = listmatlabformat( uout_nuc);
readlist_tub = listmatlabformat( uout_tub);
readlist_er = listmatlabformat( uout_er);
% $$$ if strcmpi(naming_convention.nuclear_channel, 'ignore_channel')
% $$$   readlist_nuc = repmat({''}, size(readlist)); 
% $$$ else
% $$$   uout_nuc = findreplacestring(uout, naming_convention.protein_channel, naming_convention.nuclear_channel);
% $$$   readlist_nuc = listmatlabformat( uout_nuc);
% $$$ end
% $$$ uout_tub = findreplacestring(uout, naming_convention.protein_channel, naming_convention.tubulin_channel);
% $$$ uout_er = findreplacestring(uout, naming_convention.protein_channel, naming_convention.er_channel);
% $$$ readlist_tub = listmatlabformat( uout_tub);
% $$$ readlist_er = listmatlabformat( uout_er);


mout = uout;
mout = findreplacestring( mout, naming_convention.protein_channel, naming_convention.segmentation_suffix);
%mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '/', '_');
mout = findreplacestring( mout, readdir_,writedir);
mout = findreplacestring( mout, '.tif','.png');

writelist = listmatlabformat( mout);

mkdir('./tmp')
for i=1:length(readlist)
    i
    if exist(writelist{i},'file')
        continue;
    end

    tmpfile = writelist{i};
    tmpfile(find(tmpfile=='/')) = [];
    tmpfile(find(tmpfile=='.')) = [];
    tmpfile = ['./tmp/' tmpfile '.txt'];

    if exist(tmpfile,'file')
        continue;
    end
    fid = fopen(tmpfile,'w');

    nucim = imread(readlist_nuc{i});
    cellim = imread(readlist_er{i});

    regions = segmentation( nucim, cellim, MINNUCLEUSDIAMETER, MAXNUCLEUSDIAMETER, IMAGEPIXELSIZE);

    % saving results
    imwrite( regions, writelist{i});

    fclose(fid);
    delete(tmpfile);
end
