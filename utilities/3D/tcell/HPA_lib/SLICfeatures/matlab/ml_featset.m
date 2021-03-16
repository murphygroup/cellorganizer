function [names, features, slfnames] = ml_featset( image, cropimage, ...
						  dnaimage, featsetname, ...
						  scale, radius, ...
						  bgsub, threshmeth, ...
						  har_pixsize, har_intbins)
% ML_FEATSET    Feature calculation for raw images
% [NAMES, FEATURES, SLFNAMES] = ML_FEATSET( IMAGE, CROPIMAGE, DNAIMAGE,
%                                           FEATSETNAME, SCALE, RADIUS, BGSUB)
%
% Calculates sets of features such as all 84 features, no_dna
% features 37 best features etc (see below for a full list).
% The features are returned in FEATURES and the corresponding
% feature names in NAMES. If after preprocessing all image pixels
% are zero then a vector of zeros will be returned for FEATURES and
% an empty matrix for NAMES.
% IMAGE is the raw image. This image is required.
% CROPIMAGE defines the region of interest. This is optional. If
% you don't have one, supply the empty matrix [] as an argument.
% DNAIMAGE is used to calculate DNA related features. This is
% optional. If you don't have one, supply the empty matrix []
% as an argument.
% FEATSETNAME identifies the feature set to be calculated (and num):
%             'best37'        - SLF5 (37)
%             'best34'        - Like SLF5 without DNA (34)
%             'all84'         - SLF4 (84)
%             'all114'        - SLF + wavelet features (114)
%             'allbutDNA108'  - SLF3 + wavelet features (108)
%             'nodna78'       - SLF3 (no dna) (78)
%             'nodnanohull75' - SLF4 without hull features (75)
%             'nodna_nohar'   - SLF6 (65)
%             'allwav'        - Wavelet features (30)
%             'SLF7_wav'      - SLF7 + wavelet features (114)
%             'SLF7_noskel'   - SLF7 without skeleton features (79)
%             'skel'          - Skeleton features (5)
%             'SLF7'          - SLF7 (84)
%             'SLF7dna'       - SLF7 + DNA features (90)
%             'SLF8'          - SLF8 (32)
%             'SLF8dna'       - SLF8 + DNA features (42)
%             'SLF13'         - SLF13 (31)
%             'edge'          - Edge features (5)
%             'har'           - Haralick texture features (13)
%             'gab'         - Gabor Wavelet features (60)
%             'zer'           - Zernike moment features (49)
%             'zerhar'        - zer + har (62)
%             'imgdna'        - SLF2 without hull and edge features (14)
%             'alpha'         - alpha hull features (110)
% Prepending 'mb_' to the FEATSETNAME will cause the preprocessing
% to be done the Michael Boland way. E.g. 'mb_all84' will give all
% 84 features calculated exactly the Boland way. The difference is
% in the order of thresholding/cropping.
% SCALE is micrometers per pixel. If 0 or [] is given, the default scale
% of 0.23 um/pixel is used.
% RADIUS is the cell radius in micrometers. If 0 or [] is given, the
% default radius of 34.5 um is used
% BGSUB should be omitted if you want background subtraction, or
% 'nobgsub' if you do not want it.
% THRESHMETH method of thresholding
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

% By Meel Velliste, Summer 2000
% Modified 1/16/02 MV:
%       1) Corrected error message about how many arguments are required.
%       2) Changed the help text to reflect the change in ml_features
%          which allows [] to be entered to get default value for 'scale'
%          or 'radius'.
%
% Modified 3/7/03 EJSR:
%       Formatting, comments, ml_ified
% Modified 3/9/03 EJSR:
%       SLF7_noskel index
%       best37 (idx_sasbest2) inbuilt
%       best34 (idx_34best_from78) inbuilt
%       removed loading of ml_feat_subset_idx

% Check arguments & use defaults where necessary
MaxArgs = 10;
RequiredArgs = 4;
DEFAULT_HAR_BIN = 256;
DEFAULT_HAR_PIXSIZE = 1.15;

if( nargin <= MaxArgs-1)
        har_intbins = DEFAULT_HAR_BIN;
end
if( nargin <= MaxArgs-2)
        har_pixsize = DEFAULT_HAR_PIXSIZE;
end
if( nargin <= MaxArgs-3)
    threshmeth = 'nih';
end
if( nargin <= MaxArgs-4)
        bgsub = 'yesbgsub';
end
if( nargin <= MaxArgs-5) % two missing argument
	radius = 0;
end
if( nargin <= MaxArgs-6) % three missing arguments
    scale = 0;
end
if( nargin < RequiredArgs) % four or more missing args
    error(['This function takes at least' sprintf('%i',RequiredArgs) ' arguments']);
end

if nargin<MaxArgs-3
    threshmeth = 'nih';
end

% Sort out which features need calculating
if iscell(featsetname)
    subsetlist=featsetname;
else
    % Determine if user requested Michael Boland way of preprocessing
    if( strcmp(featsetname(1:3),'mb_'))
        way = 'mb';
        featsetname(1:3) = [];
    else
        way = 'ml';
    end
  
    switch featsetname
    case 'best37'
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 1;
        subset_index = [ 3 24 14 20 78 6 2 74 15 80 11 38 61 60 79 69 34 ...
                25 5 82 9 67 68 73 72 1 76 8 48 19 7 54 4 70 ...
                77 71 42];
    case 'best34'
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [ 3 18 14 72 6 2 68 9 74 32 55 54 73 63 28 19 ...
                5 76 61 62 67 66 1 70 8 42 13 7 48 4 64 71 ...
                65 36];
    case 'all84'
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 1;
        subset_index = [1:84];
    case 'all114'
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg' 'wav'};
        use_dna = 1;
        subset_index = [1:114];
    case {'allbutDNA108','allbutdna108'}
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg' 'wav'};
        use_dna = 0;
        subset_index = [1:108];
    case {'nodna78','noDNA78'}
        subsetlist = {'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [1:78];
    case {'nodnanohull75','noDNAnohull75'}
        subsetlist = {'img' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [1:75];
    case {'nodna_nohar','noDNA_nohar'}
        subsetlist = {'img' 'hul' 'zer' 'edg'};
        use_dna = 0;
        subset_index = [1:65];
    case 'allwav'
        subsetlist = {'wav'};
        use_dna = 0;
        subset_index = [1:30];
    case 'SLF7_wav'
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg' 'wav'};
        use_dna = 0;
        subset_index = [1:114];
    case {'SLF7dna_wav','SLF7DNA_wav'}
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg' 'wav'};
        use_dna = 1;
        subset_index = [1:120];
    case 'all180'
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg' 'wav' 'gab'};
        use_dna = 1;
        subset_index = [1:180];
    case 'SLF7_noskel'
        subsetlist = {'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [1:79];
    case 'SLF7DNA_noskel'
        subsetlist = {'nof' 'img' 'hul' 'zer' 'har' 'edg'};
	use_dna = 1;
	subset_index = [1:85];
    case 'skel'
        subsetlist = {'skl'};
        use_dna = 0;
        subset_index = [1:5];
    case 'SLF7'
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [1:84];
    case {'SLF7dna','SLF7DNA'}
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 1;
        subset_index = [1:90];
    case 'SLF8'
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 0;
        subset_index = [ 9 75 20 6 72 77 24 80 8 12 69 60 14 82 48 71 3 ...
                7 25 67 1 70 51 11 5 78 81 74 27 79 73 13];
    case {'SLF8dna','SLF8DNA'}
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 1;
        subset_index = [ 9 81 20 26 6 78 30 83 86 75 12 17 88 8 44 ...
                31 3 67 19 77 7 57 84 11 73 15 5 14 87 76 ...
                74 13 40 42 28 4 79 80 54 1 49 16];
    case 'SLF13'
        subsetlist = {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg'};
        use_dna = 1;
        subset_index = [ 9 81 20 26 6 78 30 83 86 75 12 17 88 8 44 ...
                31 3 67 19 77 7 57 84 11 73 15 5 14 87 76 ...
                74];
    case 'edge'
        subsetlist = {'edg'};
        use_dna = 0;
        subset_index = [1:5];
    case 'har'
        subsetlist = {'har'};
        use_dna = 0;
        subset_index = [1:13];
    case 'gabor'
        subsetlist = {'gab'};
        use_dna = 0;
        subset_index = [1:60];
    case 'zer'
        subsetlist = {'zer'};
        use_dna = 0;
        subset_index = [1:49];
    case 'zerhar'
        subsetlist = {'zer','har'};
        use_dna = 0;
        subset_index = [1:49 50:62];
    case {'imgdna','imgDNA'}
        subsetlist = {'img'};
        use_dna = 1;
        subset_index = [1:14];
    case 'imgfield'
        subsetlist = {'img'};
	use_dna = 0;
	subset_index = [3:5];
    case 'edghar'
        subsetlist = {'edg','har'};
	use_dna = 0;
	subset_index = [1:18];
    case 'alpha'
        subsetlist = {'alpha'};
        use_dna = 0;
        subset_index = [1:110];
        % multi-cell (field level) feature set
    case 'mcell_mor'
        subsetlist = {'img','edg','skl'};
        use_dna = 0;
        %tz- 02-Apr-2006
        %subset_index=[1:18];
        %tz--
        
        %tz+ 02-Apr-2006
        subset_index=[3 4 5 9:18];
        %tz++
    case 'mcell_all'
        subsetlist = {'gab','har','img','edg','skl'};
        use_dna = 0;
        %tz- 02-Apr-2006
        %subset_index = [1:(60+13+18)];
        %tz--
        
        %tz+ 02-Apr-2006
        subset_index=[1:60 61:73 76 77 78 82:91];
        %tz++
    case 'mcell'
        subsetlist = {'har','img','edg','skl'};
        use_dna = 0;
        %modified by Juchang 04.26.06
	%subset_index = [1:(13+18)];
	subset_index = [1:13 16:18 22:31];
    case 'mcell_noskel'
        subsetlist = {'har','img','edg'};
        use_dna = 0;
        subset_index = [1:13 16:18 22:26];
    % end of modification by Juchang	
    otherwise
        error( 'Unrecognized feature set name');
    end
end

% Preprocess the image
[procimage, prot_maskimage,nonobjimg] = ml_preprocess( image, cropimage, ...
						  way, bgsub,threshmeth);
if( isempty(procimage))
  features = zeros(1,length(subset_index))+NaN;
  slfnames = [];
  names = [];
  return;
end

% Preprocess DNA image if necessary
if( use_dna)
    if( isempty( dnaimage))
        error( 'DNA image is required for this feature set.');
    end
    [procdnaimage, dna_maskimage] = ml_preprocess( dnaimage, ...
        cropimage, way,'yesbgsub',threshmeth);
    if( isempty(procdnaimage))
        features = zeros(1,length(subset_index))+NaN;
        slfnames = [];
        names = [];
        return
    end
else
    procdnaimage = [];
end

% Calculate features
[n,f,s] = ml_features( procimage, procdnaimage, prot_maskimage, ...
    subsetlist, scale, radius, nonobjimg, har_pixsize, har_intbins);

%t+
if isempty(f)
    features = zeros(1,length(subset_index))+NaN;
    names = {};
    slfnames = {};
else
% t++
    % Output features
    features = f(subset_index);
    names = n(subset_index);
    slfnames = s(subset_index);
%t+
end
%t++