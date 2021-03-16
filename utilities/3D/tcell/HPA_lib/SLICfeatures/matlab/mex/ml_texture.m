function hfeatures = ml_texture(I)
% ML_TEXTURE(I) Haralick texture features for image I
% V = ML_TEXTURE(I),
%     Returns an array of texture features for image I.
%     The value for each of the following 13 statistics is
%     included for each of the four directions.  The mean and range of
%     the features over these four directions are also included 
%     (i.e. for each feature -- 0 deg, 45 deg, 90 deg, 135 deg, mean, range)
%
%     1) Angular second moment
%     2) Contrast
%     3) Correlation
%     4) Sum of squares
%     5) Inverse difference moment
%     6) Sum average
%     7) Sum variance
%     8) Sum entropy
%     9) Entropy
%    10) Difference variance
%    11) Difference entropy
%    12) Information measure of correlation 1
%    13) Information measure of correlation 2
%    14) Maximal correlation coefficient - NOT calculated.
%
%    Reference - Haralick, RM, Shanmugam, K, Dinstein, I. (1973)  
%      Textural Features for Image Classification.  IEEE Trans.
%      on Systems, Man, and Cybernetics.  SMC-3(6):610-623.
%
%

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

% 19 Dec 98 - M.V. Boland
%    
