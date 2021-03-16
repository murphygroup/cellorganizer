function img_skel = ml_mmthin(bin_image)

%ML_MMTHIN image transformation by thinning
% IMG_SKEL=ML_MMTHIN(BIN_IMAGE)
%
% rewrite the mmthin function in the morphological toolbox
% This code is written by yenixsa and Sam in Summer 2004
% Last updated on 12/3/2005

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

degrees = 45;
num_elem = abs(360/degrees);

struct_elem{1} = [0 0 0; 2 1 2; 1 1 1];
struct_elem{2} = [2 0 0; 1 1 0; 1 1 2];
struct_elem{3} = [1 2 0; 1 1 0; 1 2 0];
struct_elem{4} = [1 1 2; 1 1 0; 2 0 0];
struct_elem{5} = [1 1 1; 2 1 2; 0 0 0];
struct_elem{6} = [2 1 1; 0 1 1; 0 0 2];
struct_elem{7} = [0 2 1; 0 1 1; 0 2 1];
struct_elem{8} = [0 0 2; 0 1 1; 2 1 1];

[r, c] = size(bin_image);

total_oper = 0;
acnum_elem = 0;
keep_num = 0;

bin_image = double(bin_image);

while 1
   acnum_elem = acnum_elem + 1;
   
   image_exp = zeros(r+4, c+4);
   image_exp(3:r+2, 3:c+2) = bin_image;
  
   [changed_image, oper] = ml_hitmiss(struct_elem{acnum_elem}, image_exp);

   total_oper = total_oper + oper;
   keep_num = acnum_elem;
    
    if(acnum_elem == num_elem)
       acnum_elem = 0;        
    end
     
    bin_image = bin_image - changed_image;
    
    if(acnum_elem == 0)
	if total_oper == 0
	    break;
        end
	total_oper = 0;
    end
end

img_skel = bin_image;
