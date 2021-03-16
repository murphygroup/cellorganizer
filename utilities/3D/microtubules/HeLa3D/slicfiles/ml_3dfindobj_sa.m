function objects = ml_3dfindobj_sa(binimg,findholes,min_obj_size,tmpdir)

% OBJECTS = ML_3DFINDOBJ_SA( BINIMG, FINDHOLES, MIN_OBJ_SIZE, TMPDIR)
%
% Find the objects in binary image BINIMG. If FINDHOLES is nonzero then
% holes are found too. Works the same as ml_3dfindobj, just uses a
% standalone program instead of a MEX file. The MEX program ran into
% trouble when allocating memory, probably because of a bug either in
% matlab or in malloc. If MIN_OBJ_SIZE is specified, then objects smaller
% than that size are ignored.

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

% Meel Velliste 5/20/02

if nargin<4
    tmpdir = tempdir;
end

if( tmpdir(end) ~= filesep) 
    tmpdir(end+1) = filesep; 
end;

if( ~exist( 'min_obj_size', 'var'))
    min_obj_size = 1;
end
    
tic
%binpath = 'bin';


bincmd=which('ml_3dfindobj_sa');
slashpos=find(bincmd==filesep);
bincmd=[bincmd(1:slashpos(end-1)) 'bin' bincmd(slashpos(end):end-2)];

fname1 = ml_fopentemp( tmpdir);
fname2 = ml_fopentemp( tmpdir);
full_imgfilename = [tmpdir fname1];
full_objfilename = [tmpdir fname2];

if strcmp(class(binimg),'double')
    binimg = uint8(floor(binimg));
end

ml_save3di( binimg, full_imgfilename);
cmdline = [bincmd ' ' full_imgfilename ' ' full_objfilename ...
          ' ' num2str(findholes) ' ' num2str(min_obj_size)];
[status, output] = system( cmdline);
if( status ~= 0)
    disp( output);
    error('There was an error in the standalone ml_3dfindobj executable');
end
ObjectFindingTime = toc
tic
% Read the object finding results
objects = readobjfile( full_objfilename);

% Remove temporary files
unix(['rm ' full_imgfilename]);
unix(['rm ' full_objfilename]);

ObjectFileReadingTime = toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function objlist = readobjfile(Filename)
% RETURN = XC_READ3DI(FILENAME)
% Read objects from Filename and convert them to the format 
% requied by ml_obj2feat.m
%
% Filename is the name of file contains this information
% Xiang Chen 5/20/02, modified by MV 5/20/02
% T. Zhao 11/17/05
%   * handle no object

fid = fopen(Filename);
if (fid == -1)
     error('Object file does not exsist');
end

objlist = {}; %*

% Read Number of Objects
objno = fread(fid, 1, 'uint32');
for m = 1 : objno
     % number of voxel in the object
     voxelno = fread(fid, 1, 'uint32');
     % number of holes in the object
     holeno = fread(fid, 1, 'uint32');
     % read in y, x, z for each voxel
     voxel = repmat(uint16(0), 3, voxelno);
     if (~feof(fid))
          voxel = fread(fid, [3 voxelno], 'uint16');
     else
          error('File corrupted');
     end

     % create the struct for this object
     temp = struct('size', voxelno, 'voxels', voxel, 'n_holes', holeno);

     % attach it to the result
     objlist{m} = temp;
end
objlist = objlist';
fclose(fid);
