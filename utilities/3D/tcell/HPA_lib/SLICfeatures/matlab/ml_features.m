function [feat_names, feat_vals, feat_slf] = ...
    ml_features(procimage, dnaimage, imagemask, featsets, scale, radius, ...
		nonobjimg, har_pixsize, har_intbins)

% [FEAT_NAMES, FEAT_VALS, FEAT_SLF] = ML_FEATURES( PROCIMAGE, DNAIMAGE,
% IMAGEMASK, FEATSETS, SCALE, RADIUS, NONOBJIMG, HAR_INTBINS, HAR_PIXSIZE)
%
% Calculates features for preprocessed image PROCIMAGE.
% SCALE is mircometers per pixel. If scale is given as 0 or [],
% it defaults to 0.23um/pixel. 
% RADIUS is cell radius in micrometers. If it is given as 0 or [],
% then it defaults to 34.5um (the radius used in Michael's thesis work). 
% Omitting the SCALE and RADIUS arguments will have the same effect as
% entering them as 0 or [].
% FEATSETS identifies the sets of features to be calculated. There
% is a three letter code for each feature set and they are calculated 
% in the order specified. e.g. {'img','zer','wav'} means calculate Image
% features, then Zernike features, and then Wavelet features.
% img - image features (8 or 14 features).



%       If DNAIMAGE is given, then DNA features are included.
% hul - hull features (3 features)
% edg - edge features (5 features)
% mor - morphological set (includes img, hul, edg, in this order).
%       If DNAIMAGE is given, then DNA features are included in the 
%       'img' set.
% zer - zernike features (49 features)
% har - haralick texture features (13 features)
% wav - wavelet features (30 features)
% skl - skeleton features (5 features)
% nof - non-object fluorescence feature(s) (currently 1 feature)
%
% PROCIMAGE must already have been preprocessing in whatever way
% desired. This function uses PROCIMAGE as is to calculate
% features.
% NONOBJIMG is the "Non-object image" i.e. the fluorescence left
% over after thresholding the background subtracted image. Supply
% this image if you want features calculated on non-object fluor-
% escence (and remember to specify 'nof' in the feature set list).
% HAR_INTBINS specifies the number of intensity bins to use for
% Haralick texture features. If omitted or specified as [] then
% the default value of 256 will be used.
% HAR_PIXSIZE specifies the pixel size (in micrometers) to use 
% for Haralick texture features. If omitted or given as [], then
% the default value of 1.15 will be used

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

% derived from mb_imgfeatures.pl by Michael Boland and
% jg_features.m by Jim Gajnak
% many bugs fixed and sections rewritten - R.F. Murphy - September-December 1999
% this version December 3, 1999
%
% image must be an image of type double (NOT scaled to 0-1)
% if desired, image should have background subtracted BEFORE calling this routine
%
% March-April, 2000 - Extended by Meel Velliste, based on rm_features
% Added input parameters for scale (micrometers per pixel)
% and cell radius (micrometers).
% Features are now scaled according to the provided scale value.
% Also changed argument 'adh' to 'mor' to reflect 'morphological features'
% instead of 'ad hoc' features. Hmm... I can't believe Michael
% actually called them 'ad hoc' features in his thesis. Seems kind
% of informal.
%
% Standardized Jul 20, 2000 by Meel Velliste. Now all preprocessing 
% has to be done before calling this function - providing maximum
% modularity. Also, changed the way that feature sets are selected, 
% and added calculation of DNA features and Wavelet features.
%
% Added SLF names March 26, 2001 by R.F. Murphy.
%
% Reordered outputs April 20th, 2001 by E. J. S. Roques
%
% MV 1/16/2002:
%      1) Changed the code to allow 'radius' and 'scale' to be specified
%         as [] (in addition to 0) to get the default value of
%         these arguments
% MV 1/17/2002:
%      1) Replaced mb_imgfeatures with ml_imgfeatures which is the
%         optimized version by Greg Porreca. Also changed the
%         haralick code to use the resampling developed by Greg.
%
%      2) Added skeleton features and non-object fluorescence
%      feature 
% MV 6/2/2002:
%      1) Changed default number of haralick intensity bins to 256
%      2) Added SLF feature names
% EJSR 3/2003
%      ml_ified
% Notes for modification from T. Zhao
%   t+: add
%   t-: commented out


DEFAULT_INTENSITY_BINS = 256;
DEFAULT_HAR_PIXEL_SIZE = 1.15;

% Check arguments & use defaults where necessary
if( nargin < 9 | isempty(har_intbins))
    har_intbins = DEFAULT_INTENSITY_BINS;
end
if( nargin < 8 | isempty(har_pixsize))
    har_pixsize = DEFAULT_HAR_PIXEL_SIZE;
end
if( nargin < 7)
    nonobjimg=[];
end
if( nargin < 6)
	radius = [];
end
if( nargin < 5)
	scale = [];
end

if isempty(radius)
    radius = 0;
end

if isempty(scale)
    scale = 0;
end


if( nargin < 4)
	error('This function takes at least 4 arguments');
end

if( isempty(procimage))
	error('You must supply a non-zero image array');
end

if( isempty(featsets))
        error('You must choose some feature sets to calculate');
end

% Make sure image is in 'double' format
procimage = double(procimage);

if (max(max(procimage))==0)
    %t+
    feat_slf = [];
    feat_names = [];
    feat_vals = [];
    warning('Image has no non-zero pixels in it');
    return;
    %t++

    %t-
    %   error('Image has no non-zero pixels in it');
    %t--
else

    % init output variables
    feat_slf = [];
    feat_names = [];
    feat_vals = [];

    % Calculate the scale factor
    default_scale = 0.23; %%% 0.23um/pixel
    if( scale == 0)
        scale = default_scale;
        scale_factor = 1;
    else
        scale_factor = scale / default_scale;
    end;

    for i = 1 : length( featsets)
        switch featsets{i}
            case 'har' % Texture (haralick) features
                % To calculate these features the image is first resized
                % to certain pixel size (1.15 um/pixel by default) and
                % the intensities rebinned (by default to 32 levels)
                har_rescale_factor = scale / har_pixsize;
                resized_img = imresize(procimage, har_rescale_factor, 'bilinear');

                % The if statment is added to avoid abnormal abortion of matlab when
                % resized_img is too small.
                % 16-JAN-2006 T. Zhao
                if length(resized_img(:))==1 | all(resized_img(:)==0)
                    tvalues = zeros(1,13);
                else
                    if size(resized_img,1)==1
                        resized_img = [resized_img;zeros(size(resized_img))];
                    end
                    if size(resized_img,2)==1
                        resized_img = [resized_img,zeros(size(resized_img))];
                    end

                    intens_binned_img = ...
                        uint8(floor(ml_imgscale(resized_img)*har_intbins));
                    tvalues = ml_texture( intens_binned_img);

                    %t+
                    tvalues = tvalues(1:13,5)';
                    %t++
                end

                tnames = {'angular_second_moment' 'contrast' 'correlation' ...
                    'sum_of_squares' 'inverse_diff_moment' 'sum_avg' 'sum_var' ...
                    'sum_entropy' 'entropy' 'diff_var' 'diff_entropy' ...
                    'info_measure_corr_1' 'info_measure_corr_2'};
                tslf = {};
                for t_no = 1:13
                    tslf{t_no}=['SLF3.' num2str(t_no+65)];
                end

                feat_slf = [feat_slf tslf];
                feat_names = [feat_names tnames];
                feat_vals = [feat_vals tvalues];
                %t-
                % 	  feat_vals = [feat_vals tvalues(1:13,5)'];
                %t--

            case 'zer' % Zernike moment features
                % if no radius was given use default
                if( radius == 0)
                    % The default is the radius used in MVB's thesis work
                    radius = 34.5; %%% (34.5um = 150 pixels)
                end
                % Convert radius from micrometers to pixels
                radius = radius / scale;

                % ml_zernike(<image>, <degree>, <radius>);
                % <image> the image you want the features from
                % <degree> Degree through which the Zernike moments are calculated (????)
                % <radius> Zernike moments are calculated using a unit circle...this defines
                %    the number of pixels to call a unit circle...
                [znames, zvalues] = ml_zernike(procimage, 12, radius);

                zslf = {};
                for z_no = 1:49
                    zslf{z_no}=['SLF3.' num2str(z_no+16)];
                end
                feat_slf = [feat_slf zslf];
                feat_names = [feat_names znames];
                feat_vals = [feat_vals abs(zvalues)];

            case 'mor' % Morphological set
                % ml_imgfeatures(<image>, <dnaproc>)
                % ml_hullfeatures(<image>, <imagehull>)
                % ml_imgedgefeatures(<image>)
                % <image> the image you want features from
                % <dnaproc> corresponding DNA image -- not appicable, hence, the forced []
                % <imagehull> the image of the convex hull of the imagemask

                % Find convex hull (this is for hull features only)
                imagehull = ml_imgconvhull(imagemask);

                [imagenames, imagevalues, imageslf] = ml_imgfeatures(procimage, dnaimage);
                [hullnames, hullvalues, hullslf] = ml_hullfeatures(procimage, imagehull);
                [edgenames, edgevalues, edgeslf] = ml_imgedgefeatures(procimage);

                % Do the rescaling of features that are affected by image
                % scale
                if( scale_factor ~= 1)
                    %%%%%%%%%%%%%%%%You must also scale some DNA features
                    index = find( strcmp( imagenames, cellstr('object_size:average')));
                    imagevalues(index) = imagevalues(index) * scale_factor^2;
                    index = find( strcmp( imagenames, cellstr('object_size:variance')));
                    imagevalues(index) = imagevalues(index) * scale_factor^4;
                    index = find( strcmp( imagenames, cellstr('object_distance:average')));
                    imagevalues(index) = imagevalues(index) * scale_factor;
                    index = find( strcmp( imagenames, cellstr('object_distance:variance')));
                    imagevalues(index) = imagevalues(index) * scale_factor^2;
                end

                feat_slf = [feat_slf imageslf hullslf edgeslf];
                feat_names = [feat_names imagenames hullnames edgenames];
                feat_vals = [feat_vals imagevalues hullvalues edgevalues];
            case 'img'
                [imagenames, imagevalues, imageslf] = ml_imgfeatures(procimage, dnaimage);
                % Do the rescaling of features that are affected by image
                % scale
                if( scale_factor ~= 1)
                    %%%%%%%%%%%%%%%%You must also scale some DNA features
                    index = find( strcmp( imagenames, cellstr('object_size:average')));
                    imagevalues(index) = imagevalues(index) * scale_factor^2;
                    index = find( strcmp( imagenames, cellstr('object_size:variance')));
                    imagevalues(index) = imagevalues(index) * scale_factor^4;
                    index = find( strcmp( imagenames, cellstr('object_distance:average')));
                    imagevalues(index) = imagevalues(index) * scale_factor;
                    index = find( strcmp( imagenames, cellstr('object_distance:variance')));
                    imagevalues(index) = imagevalues(index) * scale_factor^2;
                end

                feat_slf = [feat_slf imageslf];
                feat_names = [feat_names imagenames];
                feat_vals = [feat_vals imagevalues];
            case 'hul'
                % Find the convex hull
                imagehull = ml_imgconvhull(imagemask);
                % Calc features
                [hullnames, hullvalues, hullslf] = ml_hullfeatures(procimage, imagehull);
                feat_slf = [feat_slf hullslf];
                feat_names = [feat_names hullnames];
                feat_vals = [feat_vals hullvalues];
            case 'alpha'
                error('Sorry, alpha features are not available in this version.');
                size(imagemask)
                size(procimage)
                filename = ['/imaging/khuang/alphafeats/' imagename '.alpha']
                f = fopen(filename, 'w');
                [xi, yi] = find(imagemask);
                lenx = size(xi, 1)

                zdist = [0:15]
                for dp = 1:length(zdist)
                    for p=1:lenx
                        fprintf(f, '%i %i %s\n',xi(p),yi(p),num2str(zdist(dp)));
                    end
                end

                fclose(f);

                unix(['/imaging/velliste/alpha/alpha-4.1/bin/delcx ' filename ' > /imaging/khuang/alphafeats/delcx.out']);
                unix(['/imaging/velliste/alpha/alpha-4.1/bin/mkalf ' filename ' > /imaging/khuang/alphafeats/mkalf.out']);

                fprintf(1, 'coordinates written.\n');

                which mv_alpha
                [alphas, sigs, signames] = mv_alpha(filename);
                which mv_alphafeats
                [alfafeats, alfafeat_names] = mv_alphafeats( alphas, sigs, signames);

                feat_vals = [feat_vals alfafeats];
                feat_names = [feat_names alfafeat_names];
                [status, output] = unix(['rm -f ' filename]);
                if( status ~= 0),status,output,error('Error in shell command'),end
                [status, output] = unix(['rm -f ' filename '.dt']);
                if( status ~= 0),status,output,error('Error in shell command'),end
                [status, output] = unix(['rm -f ' filename '.info']);
                if( status ~= 0),status,output,error('Error in shell command'),end
                [status, output] = unix(['rm -f ' filename '.alf']);
                if( status ~= 0),status,output,error('Error in shell command'),end

                feat_slf = [feat_slf 1:220];

            case 'edg'

                [edgenames, edgevalues, edgeslf] = ml_imgedgefeatures(procimage);
                feat_slf = [feat_slf edgeslf];
                feat_names = [feat_names edgenames];
                feat_vals = [feat_vals edgevalues];
            case 'wav'
                % Rotational normalization
                normimage = ml_normalize( procimage);
                % Find wavelet features
                wavevals = ml_wavefeatures( normimage);
                wavenames = {'W_h1' 'W_v1' 'W_d1' ...
                    'W_h2' 'W_v2' 'W_d2' ...
                    'W_h3' 'W_v3' 'W_d3' ...
                    'W_h4' 'W_v4' 'W_d4' ...
                    'W_h5' 'W_v5' 'W_d5' ...
                    'W_h6' 'W_v6' 'W_d6' ...
                    'W_h7' 'W_v7' 'W_d7' ...
                    'W_h8' 'W_v8' 'W_d8' ...
                    'W_h9' 'W_v9' 'W_d9' ...
                    'W_h10' 'W_v10' 'W_d10' ...
                    };
                feat_slf = [feat_slf wavenames];
                feat_names = [feat_names wavenames];
                feat_vals = [feat_vals wavevals];
            case 'skl' % Skeleton features
                [skelnames, skelvals, skelslf] = ml_imgskelfeats( procimage);
                feat_slf = [feat_slf skelslf];
                feat_names = [feat_names skelnames];
                feat_vals = [feat_vals skelvals];

            case 'gab' % Gabor wavelet features

                for i=1:60
                    gabor_names(i) = {['gabor' num2str(i)]};
                end


                normimage = ml_normalize( procimage);
                tmpdir=tempdir;
                gabor_feat=ml_gaborfeat(uint8(floor(mat2gray(normimage)*255)),tmpdir(1:end-1));

                feat_vals = [feat_vals gabor_feat];

                feat_names = [feat_names gabor_names(1:60)];
                feat_slf = [feat_slf gabor_names];
                %keyboard;

            case 'nof' % Non-object fluorescence feature(s)
                if( isempty(nonobjimg))
                    error(['Need NONOBJIMG, the image representing non-object' ...
                        ' fluorescence']);
                else
                    obj_fluor = sum(procimage(:));
                    nonobj_fluor = sum(nonobjimg(:));
                    nonobj_feat = nonobj_fluor / (obj_fluor + nonobj_fluor);
                    feat_slf = [feat_slf cellstr('SLF7.79')];
                    feat_names = [feat_names cellstr('fract_nonobj_fluor')];
                    feat_vals = [feat_vals nonobj_feat];
                end
            otherwise
                error('Unrecognized feature set identifier in FEATSETS');
        end
    end

end

