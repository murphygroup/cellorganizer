function [segdna,segcell] = seg_cell_and_nucleus( dna_image_file, ...
    cell_image_file, ...
    prot_image_file, ...
    mask_image_file, ...
    downsample, ...
    display, ...
    param, ...
    currfile )

% Robert F. Murphy (murphy@cmu.edu)
%
% Copyright (C) 2012-2014 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 21, 2012 R.F. Murphy Fix downsampling; add display param;
%                            Only keep biggest object in each slice
% March 22, 2012 R.F. Murphy & Devin Sullivan Find top and bottom slice
%                            separately for cell and nucleus
% March 23, 2012 I. Cao-Berg Fixed a bug when getting cellnum from imgdir
% March 25, 2012 I. Cao-Berg The method itself checks whether a temporary file exists
% March 27, 2012 I. Cao-Berg Removed top and bottom blank slices from DNA and cell images
% April 10, 2012 I. Cao-Berg Added parameter structure to hold debug and verbose flags
% July 23, 2012  R.F. Murphy Add param to enforce assumption that slice with
%                           largest area is the bottom (appropriate for adherent cells
% August 1, 2012 I. Cao-Berg Modified the code so that is saves the bottom slice index
%                            and top slice index of the original image with
%                            respect to segcell
% August 2, 2012 D. Sullivan Fixed a bug where the DNA top/bot where saved instead of the
%                            cell membrane
% August 30, 2012 G. Johnson Added check_img_dirs to use pre-existing directories to
%                            infer values for missing directories.
%                            Added a check to skip processing for any
%                            binary images.
%                            Moved all processing steps (deconvolution, cropping, etc) to
%                            a single routine for ease of upkeep.
% Sept 7, 2012 G Johnson     Added t/f flag for check_img_dirs
% May 15, 2013 I. Cao-Berg   Updated method to support wildcards
% July 2, 2013 R.F. Murphy   Replaced calls to "boolean" with "logical"
% July 13, 2013, G. Johnson  Changed logic such that ifthe cell image is 
%                            logical, use the protein image for the nuclear
%                            segmentation
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

% temporaryFolder = [ pwd filesep 'temp' filesep 'preprocessing' ];
temporaryFolder = [ param.preprocessingFolder ];
if ~exist( temporaryFolder, 'dir' )
    mkdir( temporaryFolder );
end

segdna = [];
segcell = [];
sumslices = [];

if nargin == 6
    param = [];
    verbose = false;
    debug = false;
elseif nargin > 8
    warning('CellOrganizer: Wrong number of input arguments.');
    return
else
    try
        verbose = param.verbose;
        if ~islogical( verbose )
            verbose = false;
        end
    catch
        verbose = false;
    end
    
    try
        debug = param.debug;
        if ~islogical( debug )
            debug = false;
        end
    catch
        debug = false;
    end
    
    try
        adherent = param.adherent;
        if ~islogical( adherent )
            adherent = true;
        end
    catch
        adherent = true;
    end
    
end

%icaoberg 26/1/2014
%updated code so that it uses the display variable from parameter structure
%rather than the display variable. i kept the same method call to allow
%backward compatibility
if isfield( param, 'display' )
    display = param.display;
else
    display = false;
end

%icaoberg fixed bug
%gj aug 29, 2012
%changed img=='/' to img==filesep
if ~exist( 'currfile', 'var' )
    try
        cellnum = param.cellnum;
    catch
        cellnum = 0;
    end
else
    cellnum = currfile;
end

%icaoberg march 26, 2012
if exist( [ temporaryFolder ...
        filesep 'cell' num2str(cellnum) '.mat' ] )
    load( [ temporaryFolder ...
        filesep 'cell' num2str(cellnum) '.mat' ] );
    return
end

%icaoberg 26/1/2014
%included these strings in the parameter structure so that preprocess can
%use them
try
    param.dna_image_file = dna_image_file;
    param.cell_image_file = cell_image_file;
    param.prot_image_file = prot_image_file;
    param.prot_image_file = mask_image_file;
catch
    disp( 'Unable to save image filenames in parameter structure' );
end

%gj aug 29, 2012
% if ~check_img_dirs(imgdir,param)
%     return;
% end

% if any(round(downsample)~=downsample)
%     warning('Downsampling vector rounded to integers');
% end
% downsample=round(downsample);


%%% Process protein image first %%%
disp( 'Processing protein image' )

%icaoberg march 12, 2014
if isa( cell_image_file, 'function_handle' )
    disp( ['Reading cell image file'] );
else
    disp( ['Reading cell image file: ' cell_image_file] );
end

try
    cellim3 = ml_readimage( cell_image_file );
catch err
%     disp( [ 'Unable to load image: ' cell_image_file ] );
    getReport( err, 'extended' )
    segdna = []; segcell = []; return
end


%icaoberg march 12, 2014
if isa( mask_image_file, 'function_handle' )
    disp( ['Reading mask image file'] );
else
    disp( ['Reading mask image file: ' mask_image_file] );
end

try
    mask = ml_readimage( mask_image_file );
catch err
    disp( [ 'Unable to load image: ' mask_image_file ] );
    getReport( err, 'extended' )
    segdna = []; segcell = []; return
end

%protim3 = tz_maskimg_3d(protim3,mask);
%D. Sullivan 6/16/13 making masks optional
% if ~isempty(mask)
%     cellim3 = tz_maskimg_3d( cellim3, mask );
% end
% clear mask

zflip = false;

%D. Sullivan 6/12/13 make sure the flag calls for this
if ismember(lower(param.train.flag),{'all', 'framework', 'microtubule'})
    %grj 3/29/13, improved this block for readability
    if length(unique(cellim3)) > 2
        %D. Sullivan 6/6/13 cellfile is never used, old name.
        cellfile=['cell_image' num2str(cellnum) '.mat'];
        %D. Sullivan 3/21/13, need to pass the psf in rather than always using
        %the HeLa one.
        %     [segcell, top_slice, bot_slice, sumslices] = preprocess(cellim3, cellfile, '3DHeLa_Cell_PSF.tif', downsample, adherent, 0.98, 0.02, display);
        if isfield(param,'psf')
            psf = param.psf;
        else
            %icaoberg
            %removing use of default 3D hela psf
            disp('No PSF specified.');
            psf = [];
        end
        
        %D.Sullivan 3/21/13 added zflip flag if the cell was upsidedown.
        %This lets us know to change the nuclear orientation
        %default is 0
        
        disp( 'Preprocessing cell image' );
        %icaoberg 26/1/2014
        %documented for debugging purposes
        
        top_thresh = 0.98; %the fraction of fluorescence you wish to consider as the top slice
        bot_thresh = 0.02; %the fraction of fluorescence you wish to consider as the bottom slice
        
        %output filename for active contour segmentation debug plots
        param.output_filename = param.cell_image_file;
        [segcell, top_slice, bot_slice, sumslices,zflip] = preprocess( cellim3, ...
            cellfile, psf, downsample, adherent, top_thresh, bot_thresh, display, mask, param );
        
    else
        %D. Sullivan 6/16/13 - force segcell to be 0/1 valued boolean.
        segcell = cellim3;
        
        %murphy 7/2/2012
        segcell = logical(segcell);
        
        %D. Sullivan 6/16/13 added downsampling here.
        fprintf(1,'%s\n','Downsampling image'); tic;
        segcell = ml_downsize(segcell,downsample,'linear'); toc
        
        a = find(sum(sum(segcell,1),2));
        bot_slice = a(1);
        top_slice = a(end);
        
        sumslices = cumsum(squeeze(sum(sum(cellim3))));
        sumslices = sumslices/sumslices(end);
    end
end
%%% Process DNA image now %%%
%protim3 = ml_loadimage([imgdir '/prot/'],'tif');
%D. Sullivan 3/21/13 added check on file structure
dnaim3 = ml_readimage( dna_image_file );

if isempty(dnaim3)
    %G. Johnson 7/13/13 if the cell image is logical, use the protein image
    %for the nuclear segmentation
    %G. Johnson 8/2/13 change from logical to contains more than 2 pixel
    %values
    disp('Finding nuclear mask')
    if length(unique(cellim3(:))) > 2
        dnaim3 = findDnaMask( ml_downsize(cellim3,downsample,'linear'), segcell, param );
    else
        %icaoberg march 12, 2014
        if isa( prot_image_file, 'function_handle' )
            disp( ['Reading protein image file'] );
        else
            disp( ['Reading protein image file: ' prot_image_file] );
        end
        
        improt = ml_readimage( prot_image_file );
        dnaim3 = findDnaMask( ml_downsize(improt,downsample,'linear'), segcell, param );
    end
end

%D. Sullivan 6/16/13 - check that there is a mask
% if ~isempty(mask)
%     dnaim3 = tz_maskimg_3d( dnaim3, mask );
% end

%%%%Preprocessing starts here
if length(unique(dnaim3)) > 2
    dnafile=['dna_image' num2str(cellnum) '.mat'];
    
    disp( 'Preprocess nuclear image' );
    param.output_filename = param.dna_image_file;
    %psf_filename = '3DHeLa_DNA_PSF.tif';
    psf_filename = '';
    [segdna, ~, ~, ~] = preprocess(dnaim3, dnafile, psf_filename, downsample, false, 0.95, 0.05, display, segcell, param );
    
    %D. Sullivan 6/12/13 moved this code since order of segmentation is now
    %different
    if zflip
        segdna = flipdim(segdna,3);
    end
    
else
    
    %D. Sullivan 3/21/13 make sure segdna is the correct size
    if all(size(dnaim3)==size(cellim3))
        %no downsampling has been done yet
        segdna = ml_downsize(dnaim3,downsample,'linear');
    elseif all(size(dnaim3) == floor(size(cellim3)./downsample))
        %all downsampling has been done previously
        segdna = dnaim3;
    else
        warning(['Preprocessed nuc image size does not match cell image size. ',...
            'This could be caused by manual cropping of nuc image. ',...
            'Forcing images to the same size. This may cause unexpected results']);
        dnadownsample = size(dnaim3)./size(cellim3);
        segdna = ml_downsize(dnaim3,dnadownsample,'linear');
    end
    
    %murphy 7/2/2013
    segdna = logical(segdna);
    
end

%Ensure that there is only one object for the cell and nuclear shape
segcell = ml_findmainobj(segcell);
segdna = ml_findmainobj(segdna);


if ismember(lower(param.train.flag),{'all','framework','microtubule'})
    original_segcellsize = size(segcell);
    segcell = segcell(:,:,bot_slice:top_slice);
    segdna = segdna(:,:,bot_slice:top_slice);
    
    % make sure that no part of the nucleus is outside the cell
    segdna = and(segdna,segcell);
end

if display
    try
%         pause
        figure;
        for i=1:size(segdna,3)
            RGB(:,:,1)=uint8(segdna(:,:,i)*255);
            RGB(:,:,2)=uint8(segcell(:,:,i)*255);
            RGB(:,:,3)=uint8(zeros(size(segdna(:,:,i))));
            imshow(RGB); pause(0.1);
        end
    catch
    end
end

%icaoberg march 26, 2012
try
    disp('Saving preprocessed files to disk.');
    
    %icaoberg 1/8/2012
    if ismember(param.train.flag,{'all', 'framework','microtubule'})
        save( [ temporaryFolder ...
            filesep 'cell' num2str(cellnum) '.mat' ], 'segcell', 'segdna', 'downsample', ...
            'original_segcellsize', 'bot_slice', 'top_slice', 'sumslices', ...
            'dna_image_file', 'cell_image_file', 'prot_image_file', 'mask_image_file' );
    else
        save( [ temporaryFolder ...
            filesep 'cell' num2str(cellnum) '.mat' ], 'segdna', 'downsample','sumslices', ...
            'dna_image_file', 'cell_image_file', 'prot_image_file', 'mask_image_file' );
    end
catch
    disp( 'Unable to save preprocessed files to disk.' );
end

end

%D. Sullivan 3/21/13 Now a separate function so it can be called outside of
%this program
% %gj aug 30, 2012
% function [segimg, top_slice, bot_slice, sumslices] = preprocess(img, imgFile, psfPath, downsample, adherent, top_thresh, bot_thresh, display)
%
% % Load PSFs
% infPSF = imfinfo( psfPath);
%
% psf = zeros( infPSF(1).Height, infPSF(1).Width, length(infPSF));
%
% for I=1:length(infPSF)
%     psf(:,:,I)=imread(psfPath,I);
% end
%
%
% psf = psf.^2; % Approximate confocal PSF
%
% % downsample PDFs and images as specified by user
% if exist(imgFile,'file')
%     load(imgFile)
% else
%     fprintf(1,'%s\n','Downsampling PSF and image'); tic;
%     psf = ml_downsize(psf,downsample,'average');
%     resizeimg = ml_downsize(img,downsample,'average'); toc
%     clear dnaim3
%     fprintf(1,'%s\n','Deconvolving image'); tic;
%     [image,psf2] = deconvblind(resizeimg,psf); toc
% %    save(dnafile,'dna_image','dnaPSF2');
% end
%
%
% fprintf(1,'%s\n','Segmenting image'); tic;
% sumslices = cumsum(squeeze(sum(sum(image))));
% sumslices = sumslices/sumslices(end);
% bot_slice=find(sumslices> bot_thresh); %0.05, 0.02
% bot_slice=bot_slice(1);
% top_slice=find(sumslices< top_thresh); %0.95, 0.98
% top_slice=top_slice(end);
% bot_slice=max(bot_slice,bot_slice);
% top_slice=min(top_slice,top_slice);
%
% fprintf(1,'Total slices=%i, Bottom slice=%i, Top slice=%i\n',length(sumslices),bot_slice,top_slice);
% segimg = active3Dsegment(image, bot_slice, top_slice,display);
%
% if adherent %false, <param>
%    [bot_slice,top_slice]=fixadherent(segimg,bot_slice,top_slice);
% end
% clear image
% toc
%
%
% end