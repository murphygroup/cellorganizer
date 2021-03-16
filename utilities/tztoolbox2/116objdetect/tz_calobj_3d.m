function protobj=tz_calobj_3d( img, mask, des, dnaimg, allprotimg, ScaleFactor, ...
				slice_ref, slice_offset)
%TZ_CALOBJ_3D Obsolete.

%function protobj=tz_calobj( img, mask, des, dnaimg, allprotimg, ScaleFactor, ...
%				slice_ref, slice_offset)
%OVERVIEW:
%   Find 3d objects
%PARAMETERS:
%   img - 3d image
%   mask - a binary mask, specifying a region of interest
%   danimg - DNA image
%   allprotimage - a parallel total protein image
%   ScaleFactor - physical distance for smallest intervals 3x1
%   slice_ref
%   slice_offset
%DESCRIPTION
%   Modified from Meel's code
% 
%HISTORY:
%   ??-???-???? Initial write TINGZ

error(tz_genmsg('of','tz_calobj_3d'));

% Enforce uint8 datatype
imgclass=class(img);
switch( imgclass)
 case 'uint8',
 otherwise,
     img=double(img);
     max_pixel = max(img(:));
     img=img/(max_pixel/256);
     img=uint8(img);
     img = mv_3dbgsub(img);
     % Mask the images
%    mask_img_fullpath = [dirname2 '/' image_name '.mat']
%    load ( mask_img_fullpath);
     no_of_slices = size( img, 3);
     for slice_no = 1:no_of_slices
         slice = img(:,:,slice_no);
         slice(find(mask==0))=0;
         img(:,:,slice_no) = slice;
     end
     %otherwise, error('IMG must be of class uint8');
end
if( ~isempty(dnaimg))
    dnaclass=class(dnaimg);
    switch( dnaclass)
     case 'uint8',
     otherwise, error('DNAIMG must be of class uint8');
    end
end
if( ~isempty(allprotimg))
    allprotclass=class(allprotimg);
    switch( allprotclass)
     case 'uint8',
     otherwise, error('ALLPROTIMG must be of class uint8');
    end
end
%if( ~isempty(mask))
%    maskclass=class(mask);
%    switch( maskclass)
%     case 'uint8',
%     otherwise, error('MASK must be of class uint8');
%    end
%end

if( ~exist('slice_ref','var'))
    slice_ref = [];
end

% Background subtract and Crop the images
%img = mv_3dbgsub( img);
%protclean = mv_mask( img, mask);
protclean = img;
if( ~isempty( dnaimg))
    %dnaimg = mv_3dbgsub( dnaimg);
    %dnaclean = mv_mask( dnaimg, mask);
    dnaclean = dnaimg;
end
if( ~isempty( slice_ref))
    % assign a reference image (DNA or prot),
    % which is used to select a 2D slice
    switch( slice_ref)
     case 'Prot', refimg = protclean;
     case 'DNA', refimg = dnaclean;;
     otherwise, error('Invalid value for SLICE_REF');
    end
    refCOF = mv_findCOF( mv_sparse( refimg));
    % Find an appropriate slice of the 3D image
    slice_no = round(refCOF(3)) + slice_offset;
    if( slice_no < 1) slice_no = 1; end
    maxZ = size(refimg,3);
    if( slice_no > maxZ) slice_no = maxZ; end
    if( ~isempty( dnaimg))
	dnaclean = dnaclean(:,:,slice_no);
    end
    protclean = protclean(:,:,slice_no);
end
if( ~isempty( dnaimg))
    dnaCOF = mv_findCOF( mv_sparse( dnaclean));
    dnathresh = 255*mb_nihthreshold( dnaclean);
    dnabin = mv_binarize( dnaclean, uint8(dnathresh));
    clear dnaclean;
else
    dnaCOF = [];
    dnabin = [];
end

cellCOF=[]; cellbin=[];
%Prot
protthresh = 255*mb_nihthreshold( protclean);
protbin = mv_binarize( protclean, uint8(protthresh));
% Calculate the features
%protobj = mv_3dfindobj( mv_majfilt(protbin), 1);
ProtObjs=mv_3dfindobj_sa( protbin,1);
[ProtCOF, ProtCOFs] = mv_findCOFs( ProtObjs, protclean);
ProtCOFDists = mv_eucdist( ProtCOF, ProtCOFs, ScaleFactor);

%protobj = struct('ProtCOF',ProtCOF,'ProtCOFs',ProtCOFs,'ProtCOFDists',ProtCOFDists,'ProtObjs',ProtObjs);
protobj=[];
save([des '/ProtCOF.mat'],'ProtCOF')
save([des '/ProtCOFs.mat'],'ProtCOFs')
save([des '/ProtCOFDists.mat'],'ProtCOFDists')
save([des '/ProtObjs.mat'],'ProtObjs')

