function [names, values, slfnames] = ml_imgfeatures(imageproc, dnaproc)

% [NAMES, VALUES, SLFNAMES] = ML_IMGFEATURES(IMAGEPROC, DNAPROC)
%
%    Calculates features for IMAGEPROC
%    where IMAGEPROC contains the pre-processed fluorescence image, 
%    and DNAPROC the pre-processed DNA fluorescence image.
%    Pre-processed means that the image has been cropped and had 
%    pixels of interest selected (via a threshold, for instance).
%    Use DNAPROC=[] to exclude features based on the DNA image.  
%
%    Features calculated include:
%      - Number of objects
%      - Euler number of the image (# of objects - # of holes)
%      - Average of the object sizes
%      - Variance of the object sizes
%      - Ratio of the largest object to the smallest
%      - Average of the object distances from the COF
%      - Variance of the object distances from the COF
%      - Ratio of the largest object distance to the smallest
%      - DNA: average of the object distances to the DNA COF
%      - DNA: variance of the object distances to the DNA COF
%      - DNA: ratio of the largest object distance to the smallest
%      - DNA/Image: distance of the DNA COF to the image COF
%      - DNA/Image: ratio of the DNA image area to the image area
%      - DNA/Image: fraction of image that overlaps with DNA 
%
% 10 Aug 98 - M.V. Boland
% 11 Jul 01 - Moment calculations optimized - G. Porreca
% 08 Aug 01 - Modified to call C implementation of previous
%             optimizations - G. Porreca
% Jun 2, 2002 - M.Velliste: added SLF names

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

% $Id: ml_imgfeatures.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

%
% Initialize the variables that will contain the names and
%   values of the features.
%
names = {} ;
slfnames = {} ;
values = [] ;

%
% Features from imfeature()
%
features = imfeature(double(im2bw(imageproc)), 'EulerNumber') ;

%
% Calculate the number of objects in IMAGE
%
imagelabeled = bwlabel(im2bw(imageproc)) ;
obj_number = max(imagelabeled(:)) ;

names = [names cellstr('object:number') cellstr('object:EulerNumber')] ;
slfnames = [slfnames cellstr('SLF1.1') cellstr('SLF1.2')] ;
values = [values obj_number features.EulerNumber] ;

%
% Calculate the center of fluorescence of IMAGE
%
imageproc_m00 = ml_imgmoments(imageproc,0,0) ;
imageproc_m01 = ml_imgmoments(imageproc,0,1) ;
imageproc_m10 = ml_imgmoments(imageproc,1,0) ;
imageproc_center = [imageproc_m10/imageproc_m00 imageproc_m01/imageproc_m00] ;

% 
% Calculate DNA COF, if necessary
%
if ~isempty(dnaproc)
	dnaproc_m00 = ml_imgmoments(dnaproc,0,0) ;
	dnaproc_m01 = ml_imgmoments(dnaproc,0,1) ;
	dnaproc_m10 = ml_imgmoments(dnaproc,1,0) ;
	dnaproc_center = [dnaproc_m10/dnaproc_m00 ...
                          dnaproc_m01/dnaproc_m00] ;
end

%
% Find the maximum and minimum object sizes, and the distance 
%    of each object to the center of fluorescence
%
obj_minsize = realmax ;
obj_maxsize = 0 ;
obj_sizes = [] ;
obj_mindist = realmax ;
obj_maxdist = 0 ;
obj_distances = [] ;

obj_dnamindist = realmax ;
obj_dnamaxdist = 0 ;
obj_dnadistances = [] ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP 08/08/01
% Calls gp_imgmoments, which calls a C implementation of the code commented
% out below; the C version runs approx. 100X faster (10,000%)
%
[img_moment00, img_moment10, img_moment01, obj_size] = ml_moments_1(int32(imageproc),... 
								     int32(imagelabeled));

obj_sizes = obj_size(1,2:end);
% /GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP 07/11/01b
% Replaces code in for(i=1:obj_number) to calculate object moments with
% increased efficiency
%moment_length = max(max(imagelabeled)) + 1;
%img_moment00 = zeros(1,moment_length);
%img_moment10 = zeros(1,moment_length);
%img_moment01 = zeros(1,moment_length);
%obj_size = zeros(1,moment_length);
%
%
%for(x = 1 : size(imageproc,2))
%     for(y = 1 : size(imageproc,1))
%
%        moment_array_index = imagelabeled(y,x) + 1;
%        img_moment00(1,moment_array_index) = img_moment00(1,moment_array_index) + imageproc(y,x);
%        img_moment10(1,moment_array_index) = img_moment10(1,moment_array_index) + (x * imageproc(y,x));
%	 img_moment01(1,moment_array_index) = img_moment01(1,moment_array_index) + (y * imageproc(y,x));
%        obj_size(1,moment_array_index) = obj_size(1,moment_array_index) + 1;
%
%     end
%end
%
% /GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for (i=1:obj_number)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GP 07/11/01b
	%obj_size = size(find(imagelabeled==i),1) ;
	%if obj_size < obj_minsize
	%     obj_minsize = obj_size;
	%end
	%if obj_size > obj_maxsize
	%     obj_maxsize = obj_size;
	%end
        %
	%obj_sizes = [obj_sizes obj_size] ;
        %
        % /GP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GP 07/11/01b
        size_obj = size(obj_size);
        if(size_obj(2) >= (i+1))
	 if obj_size(1,i+1) < obj_minsize
		obj_minsize = obj_size(1,i+1) ;
	 end
	 if obj_size(1,i+1) > obj_maxsize
		obj_maxsize = obj_size(1,i+1) ;
	 end
        end
        % /GP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GP 07/11/01b
        obj_m00 = double(img_moment00(i+1));
        obj_m10 = double(img_moment10(i+1));
        obj_m01 = double(img_moment01(i+1));
        % /GP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GP 07/11/01a
        %/obj_img = roifilt2(0,imageproc,~(imagelabeled==i));
	%/obj_m00 = mb_imgmoments(obj_img,0,0) ;
	%/obj_m10 = mb_imgmoments(obj_img,1,0) ;
	%/obj_m01 = mb_imgmoments(obj_img,0,1) ;
        % /GP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GP 07/11/01a
	%obj_m00 = mb_imgmoments(roifilt2(0,imageproc,~(imagelabeled==i)),0,0) ;
	%obj_m10 = mb_imgmoments(roifilt2(0,imageproc,~(imagelabeled==i)),1,0) ;
	%obj_m01 = mb_imgmoments(roifilt2(0,imageproc,~(imagelabeled==i)),0,1) ;
        % /GP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	obj_center = [obj_m10/obj_m00 obj_m01/obj_m00]; 
	obj_distance = sqrt((obj_center - imageproc_center)...
                             *eye(2)*(obj_center - imageproc_center)') ;
	
	if obj_distance < obj_mindist
		obj_mindist = obj_distance ;
	end
	if obj_distance > obj_maxdist
		obj_maxdist = obj_distance ;
	end

	obj_distances = [obj_distances obj_distance] ;

	if ~isempty(dnaproc)
		obj_dnadistance = sqrt((obj_center - dnaproc_center) ...
                                      *eye(2)*(obj_center - dnaproc_center)') ;
		if obj_dnadistance < obj_dnamindist
			obj_dnamindist = obj_dnadistance ;
		end
		if obj_dnadistance > obj_dnamaxdist
			obj_dnamaxdist = obj_dnadistance ;
		end

		obj_dnadistances = [obj_dnadistances obj_dnadistance] ;
	end
end

obj_size_avg = mean(obj_sizes,2) ;
obj_size_var = var(obj_sizes) ;
obj_size_ratio = obj_maxsize/obj_minsize ;

names = [names cellstr('object_size:average') ...
		cellstr('object_size:variance') ...
		cellstr('object_size:ratio')] ;
slfnames = [slfnames cellstr('SLF1.3') ...
	    cellstr('SLF1.4') ...
	    cellstr('SLF1.5')] ;
values = [values obj_size_avg obj_size_var obj_size_ratio] ;

obj_dist_avg = mean(obj_distances) ;
obj_dist_var = var(obj_distances) ;
if obj_mindist ~= 0 
	obj_dist_ratio = obj_maxdist/obj_mindist ;
else
	obj_dist_ratio = 0 ;
end
names = [names cellstr('object_distance:average') ... 
		cellstr('object_distance:variance') ...
		cellstr('object_distance:ratio')] ;
slfnames = [slfnames cellstr('SLF1.6') ... 
		cellstr('SLF1.7') ...
		cellstr('SLF1.8')] ;
values = [values obj_dist_avg obj_dist_var obj_dist_ratio] ;

if ~isempty(dnaproc)
	obj_dnadist_avg = mean(obj_dnadistances) ;
	obj_dnadist_var = var(obj_dnadistances) ;
	if obj_dnamindist ~= 0 
		obj_dnadist_ratio = obj_dnamaxdist/obj_dnamindist ;
	else
		obj_dnadist_ratio = 0 ;
	end
	names = [names cellstr('DNA_object_distance:average') ...
                       cellstr('DNA_object_distance:variance') ...
                       cellstr('DNA_object_distance:ratio')] ;
	slfnames = [slfnames cellstr('SLF2.17') ...
                       cellstr('SLF2.18') ...
                       cellstr('SLF2.19')] ;
	values = [values obj_dnadist_avg obj_dnadist_var obj_dnadist_ratio] ;

	dna_image_distance = sqrt((imageproc_center-dnaproc_center)... 
                               *eye(2)*(imageproc_center-dnaproc_center)') ;

	dna_area = size(find(im2bw(dnaproc)),1) ;
	image_area = size(find(im2bw(imageproc)),1) ;
	%
	% what fraction of the image fluorescence area overlaps the dna image?
	%
	image_overlap = size(find(roifilt2(0,imageproc,~im2bw(dnaproc))),1) ;

	if image_area == 0
		dna_image_area_ratio = 0 ;
		image_dna_overlap = 0 ;
	else
		dna_image_area_ratio = dna_area/image_area ;
		image_dna_overlap = image_overlap/image_area ;
	end
	
	names = [names cellstr('DNA/image:distance') ...
			cellstr('DNA/image:area_ratio') ...
			cellstr('DNA/image:overlap')] ;
	slfnames = [slfnames cellstr('SLF2.20') ...
			cellstr('SLF2.21') ...
			cellstr('SLF2.22')] ;
	values = [values dna_image_distance dna_image_area_ratio ...
                  image_dna_overlap] ;
end
