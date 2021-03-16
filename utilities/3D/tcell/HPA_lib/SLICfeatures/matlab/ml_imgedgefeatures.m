function [names, values, slfnames] = ml_imgedgefeatures(imageproc)
% MB_IMGFEATURES - calculates features for IMAGEPROC
% [NAMES, VALUES, SLFNAMES] = MB_IMGFEATURES(IMAGEPROC),
%    where IMAGEPROC contains the pre-processed fluorescence image.
%    Pre-processed means that the image has been cropped and had
%    pixels of interest selected (via a threshold, for instance).
%
%   Features calculated include:
%   (The following feature descriptions were added here by T. Zhao 
%   according to the reference 
%   "R. F. Murphy, M. Velliste, and G. Porreca (2003) Robust Numerical 
%   Features for Description and Classification of Subcellular Location 
%   Patterns in Fluorescence Microscope Images. J. VLSI Sig. Proc. 35: 
%   311-321.")
%   1. Fraction of above-threshold pixels along edge
%   2. Measure of edge gradient intensity homogenerity
%   3. Measure of edge direction homogenerity 1
%   4. Measure of edge direction homogenerity 2
%   5. Measure of edge direction difference
%
% 07 Mar 99 - M.V. Boland
% 
% Corrections by M. Velliste 1/20/02
% 1) corrected the way non-zero edge pixels are selected (see below).
% 2) corrected the homogeneity feature so that it is based on a
%    histogram of edge magnitudes (see below).
% 3) Made some comments below about differences between Matlab 5 &
%    6, and provided a fix for version 5 so that it would give
%    consistent results.
%
% M.Velliste June 2, 2002: added SLF names
   
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

    
% $Id: ml_imgedgefeatures.m,v 1.4 2006/06/27 13:33:47 tingz Exp $

%
% Initialize the variables that will contain the names and
%   values of the features.
%
names = {} ;
slfnames = {} ;
values = [] ;

%
% Total area of the image that is edges
%
A = bwarea(edge(imageproc,'canny',[]))/bwarea(im2bw(imageproc)) ;

%
% Directional edge filters
%
N = [1 1 1 ; 0 0 0 ; -1 -1 -1] ;
W = [1 0 -1 ; 1 0 -1 ; 1 0 -1] ;

%
% Calculation of the gradient from two orthogonal directions
%
% Comments on differences of filter2 implementations under Matlab 5
% and 6: Under Matlab 5 sometimes filter2 can return negative
% zeros (i.e. in terms of binary representation, numbers where only
% the most significant bit is set, which means a zero with the sign
% bit set). Under Matlab 6 this has been fixed and filter2 never 
% returns negative zeros, always positive. Under matlab5 this
% caused inconsistent output below from atan2. Note that both in
% Matlab 5 and 6 atan2 behaves as follows: atan2(0,-1)=pi and
% atan2(-0,-1)=-pi. So under matlab 6 everything is OK because
% negative zeros do not occur in the first place. I.e. Matlab 6 has
% the correct implementation of filter2 and matlab 5 has an
% incorrect implementation. This does mean that the features
% calculated under version 6 are different from version 5, but they
% are correct. (Note, filter2 itself was not changed. One of the
% functions it calls was changed, probably conv2 which is a
% built-in function).
iprocN = filter2(N,imageproc) ;
iprocW = filter2(W,imageproc) ;

%
% Calculate the magnitude and direction of the gradient
%
iprocmag = sqrt(iprocN.^2 + iprocW.^2) ;
iproctheta = atan2(iprocN, iprocW) ;
% As mentioned above, under Matlab 5 the filter2 function can
% sometimes return negative zeros, and as it happens
% atan2(-0,-1)=-pi. Of course 0 and -0 are the same thing, and so
% are pi and -pi. So to get consistent results in Matlab 5, we
% change any -pi to pi. So this way the output is always the same
% for west facing edges.
iproctheta(find(iproctheta==-pi))=pi;

% Change by MV:
% Identify pixels in iprocmag that are not 0
%  (i.e. iprocN and iprocW were both 0). Before
% this was incorrectly based on identifying non-zero
% pixels in iproctheta, which does remove zero-magnitude
% edges, but it also removes edges that face exactly east.
%
nonzero = find(iprocmag);
v = iproctheta(nonzero);
v_mag = iprocmag(nonzero);

%
% Histogram the gradient directions
%
h = hist(v,8) ;

%
% max/min ratio
%
[hmax maxidx] = max(h) ;
hmin = min(h) ;

if (hmin ~= 0),
  maxminratio = hmax/hmin ;
else
  maxminratio = 0 ;
end

htmp=h ;
htmp(maxidx) = 0 ;
hnextmax = max(htmp) ;
maxnextmaxratio=hmax/hnextmax ;

%
% Difference between bins of histogram at angle and angle+pi
%  In general, objects have an equal number of pixels at an angle
%  and that angle+pi. The differences are normalized to the sum of 
%  the two directions.
%
sumdiff=sum(abs(h(1:4)-h(5:8))./abs(h(1:4)+h(5:8))) ;

%
% Measure of edge homogeneity - what fraction of edge pixels are in
%  the first two bins of the histogram.
% Change by MV: Made it be based on edge magnitude histogram. Was
% incorrectly based on edge direction histogram before.
h_mag = hist(v_mag,4);
homogeneity = sum(h_mag(1))/sum(h_mag(:)) ;

names = [names cellstr('edges:area_fraction') ...
         cellstr('edges:homogeneity') ...
	 cellstr('edges:direction_maxmin_ratio') ...
         cellstr('edges:direction_maxnextmax_ratio') ...
         cellstr('edges:direction_difference')] ;
slfnames = [slfnames cellstr('SLF1.9') ...
         cellstr('SLF1.10') ...
	 cellstr('SLF1.11') ...
         cellstr('SLF1.12') ...
         cellstr('SLF1.13')] ;
values = [A homogeneity maxminratio maxnextmaxratio sumdiff] ;
