function objects=tz_find3dobj(img,cropimg,findholes,min_obj_size,tmpdir)
%TZ_FIND3DOBJ Obsolete.

%function objs=tz_find3dobj(img,cropimg,findholes,min_obj_size,tmpdir)
error(tz_genmsg('of','tz_find3dobj'));

protimg = mv_3dbgsub(img);
% protimg = protimg( crop.croprgn.y, crop.croprgn.x, crop.croprgn.z);
protclean = tz_maskimg_3d( protimg, cropimg);
clear protimg;
protthresh = 255*mb_nihthreshold( protclean);
protbin = mv_binarize( protclean, uint8(protthresh));

if( ~exist( 'min_obj_size', 'var'))
    min_obj_size = 1;
end
    
tic
binpath = '/home/tingz/matlab/function/shared/bin';

if(~exist('tmpdir','var'))
    tempdir = '~/tmp';
end

fname1 = mv_fopentemp( tempdir);
fname2 = mv_fopentemp( tempdir);
full_imgfilename = [tempdir '/' fname1];
full_objfilename = [tempdir '/' fname2];
ml_save3di( uint8(protbin), full_imgfilename);
cmdline = [binpath '/mv_3dfindobj ' full_imgfilename ' ' full_objfilename ...
          ' ' num2str(findholes) ' ' num2str(min_obj_size)];
[status, output] = unix( cmdline);
if( status ~= 0)
    disp( output);
    error('There was an error in the standalone mv_3dfindobj executable');
end
ObjectFindingTime = toc
tic
% Read the object finding results
objects = readobjfile( full_objfilename);

for i=1:length(objects)
    objects{i}.gray=protclean(sub2ind(size(protclean),objects{i}.voxels(1,:),objects{i}.voxels(2,:),objects{i}.voxels(3,:)));
    objects{i}.imgsize=size(img);
    if any(objects{i}.gray==0)
        warning('object finding might be wrong')
        keyboard
    end
end

% Remove temporary files
unix(['rm ' full_imgfilename]);
unix(['rm ' full_objfilename]);

ObjectFileReadingTime = toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function objlist = readobjfile(Filename)
% RETURN = XC_READ3DI(FILENAME)
% Read objects from Filename and convert them to the format 
% requied by mv_obj2feat.m
%
% Filename is the name of file contains this information
% Xiang Chen 5/20/02, modified by MV 5/20/02


fid = fopen(Filename);
if (fid == -1)
     error('Object file does not exsist');
end

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
