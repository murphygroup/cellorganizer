function calcRegionFeat(rootdir, maskdir, writedir, naming_convention, fsetnames, optimize)
%function calcRegionFeat(rootdir, maskdir, writedir, naming_convention, fsetnames, optimize)
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

% 10 Jan 10 - Jieyue Li
% 2011-11-09 tebuck: adding option to use IFconfocal_CellCycle-like
% image names.
% 2011-11-15 tebuck: channel_as_protein sets the protein image to
% that channel's image and blanks the channel's image. Seems like
% loadProteinChannel is never called, so important change is
% "switch lower(channel_as_protein)" right before
% commonScriptCalculateSet.
% 2011-11-18 tebuck: made fsetnames an arg. INCOMPLETE, NEEDS
% ACCOMODATIONS LIKE WITH numfeats.
% 2011-12-07 tebuck: Adding option
% naming_convention.segmentation_suffix to allow reordering
% channels without erroneously affecting segmentation. Adding
% naming_convention.blank_channels to blank any channels for
% feature computation, not just as with channel_as_protein.
% 2011-12-28 tebuck: channel_as_protein is redundant, removing.
% 2011-12-30 tebuck: saving computation_time per feature file.


original_fsetnames = {...
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

if ~exist('fsetnames','var') || isempty(fsetnames)
  fsetnames = original_fsetnames; 
end


datasettype = 'Region';
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
if ~isfield(naming_convention, 'blank_channels')
  naming_convention.blank_channels = {}; 
end

% Optionally mark some channels as blank:
nuclear_channel_blank = false;
tubulin_channel_blank = false;
er_channel_blank = false;
protein_channel_blank = false;
for blank_channel_index = 1:length(naming_convention.blank_channels)
  switch lower(naming_convention.blank_channels{blank_channel_index})
   case 'nuclear'
    nuclear_channel_blank = true;
   case 'tubulin'
    tubulin_channel_blank = true;
   case 'er'
    er_channel_blank = true;
   case 'protein'
    protein_channel_blank = true;
   otherwise
    error('blank_channels must be a cell array of zero or more of the strings ''nuclear'', ''tubulin'', ''er'', or ''protein''')
  end
end



% $$$ if ~exist('channel_as_protein','var')
% $$$     channel_as_protein = 'protein';
% $$$ end


if ~exist('optimize','var')
    optimize = false;
end



filetype = 'tif';
%greparg = '| grep green';
greparg = ['| grep ', naming_convention.protein_channel];

ind = find(rootdir=='/');
rootdir_ = rootdir;
rootdir_(ind) = '_';

uout = unixfind( rootdir, filetype, greparg);
readlist = listmatlabformat( uout);

uout_nuc = findreplacestring(uout, naming_convention.protein_channel, naming_convention.nuclear_channel);
uout_tub = findreplacestring(uout, naming_convention.protein_channel, naming_convention.tubulin_channel);
uout_er = findreplacestring(uout, naming_convention.protein_channel, naming_convention.er_channel);
readlist_nuc = listmatlabformat( uout_nuc);
readlist_tub = listmatlabformat( uout_tub);
readlist_er = listmatlabformat( uout_er);



mout = uout;
mout = findreplacestring( mout, naming_convention.protein_channel, naming_convention.segmentation_suffix);
%mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '/', '_');
mout = findreplacestring( mout, '.tif', '.png');
%mout = findreplacestring( mout, rootdir_, ['./data/masks/']);
mout = findreplacestring( mout, rootdir_, maskdir);
readlist_mask = listmatlabformat( mout);


cleanobject.channel_path = [];
cleanobject.channel = [];
cleanobject.channel_mcp = [];
cleanobject.channel_filtered = [];
cleanobject.channel_thr = [];
cleanobject.channel_mthr = [];
cleanobject.channel_fg = [];
cleanobject.channel_regions = [];
cleanobject.channel_objectsizes = [];
cleanobject.downsampled2x = [];
cleanobject.downsampled2x_filtered = [];
cleanobject.downsampled2x_thr = [];
cleanobject.downsampled2x_mthr = [];
cleanobject.downsampled2x_fg = [];
cleanobject.downsampled2x_objectsizes = [];
cleanobject.downsampled4x = [];
DSF = 2;
mkdir('./tmp')
for i=1:length(readlist)
    i
    % 2011-12-31 tebuck: commented because it would interfere with
    % paraellel computation of multiple calls with the same
    % temporary directory:
% $$$     % Nevermind, just assume, e.g.,
% $$$     % get_concatenated_region_features is only called once per
% $$$     % directory of images (with many directories with few images
% $$$     % per directory, this should be a safe way to parallelize?):
% $$$     tmpfile = [readlist{i} '_' datasettype];
% $$$     tmpfile(find(tmpfile=='/')) = [];
% $$$     tmpfile(find(tmpfile=='.')) = [];
% $$$     tmpfile = ['./tmp/' tmpfile '_REGION.txt'];
% $$$ 
% $$$     if exist(tmpfile,'file')
% $$$         continue;
% $$$     end
% $$$     fid = fopen(tmpfile,'w');

    protfieldstruct = cleanobject;
    nucfieldstruct = cleanobject;
    tubfieldstruct = cleanobject;
    erfieldstruct = cleanobject;
    maskfieldstruct = cleanobject;

    for zed = 1:length(fsetnames)
% $$$         mout = findreplacestring( readlist{i}, '/', '_');
        mout = readlist{i};
        mout = findreplacestring( mout, naming_convention.protein_channel, naming_convention.segmentation_suffix);
        mout = findreplacestring( mout, '/', '_');
        feature_set_suffix = ['_', fsetnames{zed}];
% $$$         if ~strcmpi(channel_as_protein, 'protein')
% $$$           feature_set_suffix = [feature_set_suffix, '_', channel_as_protein, '-focus'];
% $$$         end
        mout = findreplacestring( mout, '.tif', [feature_set_suffix, '.mat']);
        writepath = findreplacestring( mout, rootdir_, [writedir 'features/' lower(datasettype) '/']); %%
        
        if exist(writepath,'file')
            continue;
        end
        
        tmpfile2 = writepath;
        tmpfile2(find(tmpfile2=='/')) = [];
        tmpfile2(find(tmpfile2=='.')) = [];
        tmpfile2 = ['./tmp/' tmpfile2 '_REGION.txt']; %%

        if exist(tmpfile2,'file')
            continue;
        end
        fid2 = fopen(tmpfile2,'w');
        
        start_time = tic;

        protfieldstruct.channel_path = readlist{i};
        nucfieldstruct.channel_path = readlist_nuc{i};
        tubfieldstruct.channel_path = readlist_tub{i};
        erfieldstruct.channel_path = readlist_er{i};
        maskfieldstruct.channel_path = readlist_mask{i};

        maskAllChannels
        
        allfeats = [];
        for j=1:length(protfieldstruct.channel_regions)
            protstruct = cleanobject;
            nucstruct = cleanobject;
            tubstruct = cleanobject;
            erstruct = cleanobject;

            protstruct.channel_path = readlist{i};
            nucstruct.channel_path = readlist_nuc{i};
            tubstruct.channel_path = readlist_tub{i};
            erstruct.channel_path = readlist_er{i};
            maskstruct.channel_path = readlist_mask{i};

            protstruct.channel = protfieldstruct.channel_regions{j};
            nucstruct.channel = nucfieldstruct.channel_regions{j};
            tubstruct.channel = tubfieldstruct.channel_regions{j};
            erstruct.channel = erfieldstruct.channel_regions{j};

            if nuclear_channel_blank
              nucstruct.channel = nucstruct.channel * 0; 
            end
            if tubulin_channel_blank
              tubstruct.channel = tubstruct.channel * 0; 
            end
            if er_channel_blank
              erstruct.channel = erstruct.channel * 0; 
            end
            if protein_channel_blank
              protstruct.channel = protstruct.channel * 0; 
            end

            commonScriptCalculateSet

            allfeats = [allfeats; feats];
        end

        feats = allfeats;
        impath = readlist{i};

        % saving results
        computation_time = toc(start_time);
        save( writepath,'feats','names','slfnames','impath', 'computation_time');

        fclose(fid2);
        delete(tmpfile2);
    end
% $$$     fclose(fid);
% $$$     delete(tmpfile);
end
