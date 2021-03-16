function spfeat = gmm_allnucspfeat( method, dnaimgdir, cellimgdir, protimgdir,  param )

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% July 23, 2012 R.F. Murphy Add param and make default verbose
% August 2, 2012 I. Cao-Berg Fixed bug where the image directory loop
%                   would ignore the first two images
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
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

if nargin == 2
   param = [];
end

try
  verbose = param.verbose;
catch
  param.verbose = true;
end

try
  debug = param.debug;
catch
  param.debug = true;
end

dna_image_files = ml_ls( dnaimgdir );
cell_image_files = ml_ls( cellimgdir );
prot_image_files = ml_ls( protimgdir );

try
    masks = ml_ls( param.masks );
catch
    masks = [];
    
end

%icaoberg 5/15/2013
spfeat = {};
for i = 1:1:length( cell_image_files )
    %icaoberg 15/5/2013
    try
        dna_image_file = dna_image_files{i};
    catch
        dna_image_file = '';
    end
    
    cell_image_file = cell_image_files{i};
    
    try
        prot_image_file = prot_image_files{i};
    catch
        prot_image_file = '';
    end
    
    try
        mask_image_file = masks{i};
    catch
        mask_image_file = [];
    end
    
    disp( ['Image:' num2str(i) '/' num2str( length( cell_image_files ) ) ] );
    if exist( dna_image_file ) || exist( cell_image_file )
        try
            param.cellnum = i;
            temp = tp_nucimgfeat( dna_image_file, ...
                cell_image_file, ...
                prot_image_file, ...
                mask_image_file, ...
                method, param );
            
            temp.dna_image_file = dna_image_file;
            temp.cell_image_file = cell_image_file;
            temp.prot_image_file = prot_image_file;
            temp.mask_image_file = mask_image_file;
            
            if ~isempty( temp )
                spfeat{length(spfeat)+1} = temp;
                disp( 'Features successfully calculated' );
            else
                continue
            end
        catch err
            % need to handle this error
            disp( 'Unable to calculate features' );
            continue
        end
    end
end
end
