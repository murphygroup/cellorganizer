function [segdna,segcell, bot_slice,top_slice] = seg_cell_and_nucleus_3D( imdna, ...
    imcell, ...
    improt, ...
    immask, ...
    downsample, ...
    options )

% Robert F. Murphy (murphy@cmu.edu)
%
% Copyright (C) 2012-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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
% March 12, 2017 R.F. Murphy Display central slice of final result

if ~exist('options', 'var')
    options = [];
end

options = ml_initparam(options, struct('train', [], ...
                                    'verbose', false, ...
                                    'debug', false, ...
                                    'adherent', true, ...
                                    'display', false));

options.train = ml_initparam(options.train, struct('flag', 'all'));

segdna = [];
segcell = [];
sumslices = [];

zflip = false;

if length(unique(imcell)) > 2
    if isfield(options,'psf')
        psf = options.psf;
    else
        %icaoberg
        %removing use of default 3D hela psf
        disp('No PSF specified.');
        psf = [];
    end

    %D.Sullivan 3/21/13 added zflip flag if the cell was upsidedown.
    %This lets us know to change the nuclear orientation
    %default is 0

    if options.verbose, disp( 'Preprocessing cell image' ), end;
    %icaoberg 26/1/2014
    %documented for debugging purposes

    top_thresh = 0.98; %the fraction of fluorescence you wish to consider as the top slice
    bot_thresh = 0.02; %the fraction of fluorescence you wish to consider as the bottom slice

    %output filename for active contour segmentation debug plots
    [segcell, top_slice, bot_slice, sumslices,zflip] = preprocess( imcell, ...
        [], psf, downsample, options.adherent, top_thresh, bot_thresh, options.display, immask, options );
    
else
    %D. Sullivan 6/16/13 - force segcell to be 0/1 valued boolean.
    segcell = imcell;

    %murphy 7/2/2012
    segcell = logical(segcell);

    %D. Sullivan 6/16/13 added downsampling here.
    if  ~all(downsample == 1)
        if options.verbose, fprintf(1,'%s\n','Downsampling image'), end;
    
        segcell = ml_downsize(segcell,downsample,'linear');
    end

    a = find(sum(sum(segcell,1),2));
    %Serena 10/20
    if isempty(a)
        bot_slice=0;
    else
        bot_slice=a(1);
    end
    top_slice = a(end);

    sumslices = cumsum(squeeze(sum(sum(imcell))));
    sumslices = sumslices/sumslices(end);
end

%%% Process DNA image now %%%
%protim3 = ml_loadimage([imgdir '/prot/'],'tif');
%D. Sullivan 3/21/13 added check on file structure

if isempty(imdna)
    %G. Johnson 7/13/13 if the cell image is logical, use the protein image
    %for the nuclear segmentation
    %G. Johnson 8/2/13 change from logical to contains more than 2 pixel
    %values
    if options.verbose, disp('Finding nuclear mask'), end;
    
    if length(unique(imcell(:))) > 2
        imdna = findDnaMask( ml_downsize(imcell,downsample,'linear'), segcell, options );
    else
        imdna = findDnaMask( ml_downsize(improt,downsample,'linear'), segcell, options );
    end
end

%D. Sullivan 6/16/13 - check that there is a mask
% if ~isempty(mask)
%     dnaim3 = tz_maskimg_3d( dnaim3, mask );
% end

%%%%Preprocessing starts here
if length(unique(imdna)) > 2
%     dnafile=['dna_image' num2str(cellnum) '.mat'];
    
    if options.verbose, disp( 'Preprocess nuclear image' ), end;
%     options.output_filename = options.dna_image_file;
    %psf_filename = '3DHeLa_DNA_PSF.tif';
    psf_filename = '';
    [segdna, ~, ~, ~] = preprocess(imdna, [], psf_filename, downsample, false, 0.95, 0.05, options.display, segcell, options );
    
    %D. Sullivan 6/12/13 moved this code since order of segmentation is now
    %different
    if zflip
        segdna = flipdim(segdna,3);
    end
    
else
    
    %D. Sullivan 3/21/13 make sure segdna is the correct size
    if all(size(imdna)==size(imcell))
        %no downsampling has been done yet
        segdna = ml_downsize(imdna,downsample,'linear');
    elseif all(size(imdna) == floor(size(imcell)./downsample))
        %all downsampling has been done previously
        segdna = imdna;
    else
        warning(['Preprocessed nuc image size does not match cell image size. ',...
            'This could be caused by manual cropping of nuc image. ',...
            'Forcing images to the same size. This may cause unexpected results']);
        dnadownsample = size(imdna)./size(imcell);
        segdna = ml_downsize(imdna,dnadownsample,'linear');
    end
    
    %murphy 7/2/2013
    segdna = logical(segdna);
    
end

%Ensure that there is only one object for the cell and nuclear shape
segcell = ml_findmainobj(segcell);
segdna = ml_findmainobj(segdna);

if ~isempty(segcell)
    segdna = and(segdna, segcell);
end

if options.display
    rm_showCentralSlice(segcell,zeros(size(segcell)),segdna);
    pause(10);
    [tpath,tname,text] = fileparts(options.cell_image_path);
    saveas(gcf, [tname 'segcellnuc'], 'fig');
end
end

function rm_showCentralSlice(r,g,b)

RGB(:,:,1) = double(r(:,:,round(end/2)))*255.;
RGB(:,:,2) = double(g(:,:,round(end/2)))*255.;
RGB(:,:,3) = double(b(:,:,round(end/2)))*255.;
imshow(RGB)
end