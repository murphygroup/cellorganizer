function passed = check_img_dirs( imgdir,param )
%CHECK_IMG_DIRS Routine to determine if any folders are missing and if so, attempt to
%make those folders and approximate the data that belongs in them (i.e. if
%the DNA folder is missing, use the cell or protein folders to find the
%"nuclear hole" and build a mask for the cell
%
% If cell dir does not exist
%   If protein directory exists
%       Use protein to find cell mask
%   Else
%       return false
%
%
% If DNA dir does not exist
%   %Use the cell image and or protein image to find the nucleus
%
                                                                   
%gj aug 29, 2012
%gj october, 2012
%gj 3/12/13 - implemented mask
%
% Copyright (C) 2013 Murphy Lab
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



err = MException('CellOrganizer:MissingImgDirs', 'Please check if all image directories exist');


%D. Sullivan 3/19/13
%can't have hard coded directories here
%Making param structure for checking the directories
if isfield(param.nucleus,'subdir')
    dnadir = [imgdir filesep param.nucleus.subdir];
else 
    warning('No dna directory specified, attempting holefinding');
    mkdir(dnadir)
end
if isfield(param.cell,'subdir')
    celldir = [imgdir filesep param.cell.subdir];
else
    error('No cell directory specified');
end

if isfield(param.protein,'subdir')
    protdir = [imgdir filesep param.protein.subdir];
else
    error('No protein directory specified');
end

if isfield(param,'mask')
    if isfield(param.mask,'subdir')
        maskdir = [imgdir filesep param.mask.subdir];
    elseif verbose==1
        warning('No optional mask directory specified');
    end
else
    maskdir = [];
end


% if isfield(param.
% protdir = [imgdir filesep 'prot' filesep];
% dnadir = [imgdir filesep 'dna' filesep];
% celldir = [imgdir filesep 'cell' filesep];
% maskdir = [imgdir filesep 'crop' filesep];

passed = true;
%if the cell dir doesnt exist
if ~exist(celldir, 'dir')
    %Use the protdir to find the cell shape
    if exist(protdir, 'dir')
        improt = ml_loadimage(protdir, 'tif');
        
        if exist(maskdir, 'dir')
            immask =  logical(ml_loadimage(maskdir, 'tif'));
            
            if size(immask,3) == 1
                immask = repmat(immask, [1,1,size(improt,3)]);
            end
        else
            immask = logical(ones(size(improt)));
        end
        
        improt(~immask) = 0;
        
        imcell = findCellMask(improt);
        
        ml_saveimage(imcell, celldir, 'tif');
    else
        passed = false;
    end
end

if ~exist(dnadir, 'dir')|| param.nucholeflag
    
    imcell = ml_loadimage(celldir, 'tif');
    
    if exist(protdir, 'dir')
        improt = ml_loadimage(protdir, 'tif');
    else
        throw(err);
    end
    
    %If the cell image is binary valued
    if length(unique(imcell)) <= 2  
        imdna = findDnaMask(improt, imcell,param);
    elseif exist(maskdir, 'dir')
        immask = ml_loadimage(maskdir, 'tif');
        imdna = findDnaMask(improt, immask,param);
    else
        imdna = findDnaMask(improt,[],param);
%         imdna = findDnaMask(imcell);
          
    end
        
    ml_saveimage(imdna, dnadir, 'tif')

end

end

