function [ seg_dna, seg_cell ] = seg_cell_and_nucleus_2D(imdna, imcell, immask, options)

% Gregory Johnson
%
% Copyright (C) 2016-2019 Murphy Lab
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

options = ml_initparam(options, ...
    struct('display', 0 ...
    ,'threshmeth_dna', 'nih' ...
    ,'threshmeth_prot', 'nih' ...
    ,'ml_objgaussmix', [] ...
    , 'dna_segment_method', 'guassian_smooth_thresh'));

if isempty(immask)
    disp('Image mask is empty')
    immask = uint8(ones(size(imdna)));
end

%  xruan 01/06/2016
% change segmentation method for imdna to make it smoother and better.
isbinary = @(x)(strcmpi(class(x),'uint8') && numel(unique(x)) == 2);

if ~isempty(imdna)
    disp('Preprocessing nuclear membrane image')
    switch options.dna_segment_method
        case 'guassian_smooth_thresh'
            if isbinary(imdna) && options.skip_preprocessing
                disp('Binary image found. Skipping preprocessing.')
                seg_dna = imdna;
            else
                procimg = ml_preprocess(double(imdna),immask,'ml','yesbgsub', options.threshmeth_dna);
                gray_img = mat2gray(double(imdna));
                gray_img = imfilter(gray_img, fspecial('gaussian', 13, 2));
                grey_thresh = graythresh(gray_img .* double(immask));
                
                if parallel.gpu.GPUDevice.isAvailable
                    gray_img = gpuArray(gray_img);
                    immask2 = gpuArray(double(immask));
                    procimg = gather(gray_img .* immask2) > grey_thresh;
                else
                    procimg = (gray_img .* double(immask)) > grey_thresh;
                end
                
                seg_dna = imfill(procimg>0,'hole');
            end
        otherwise
            if isbinary(imdna) && options.skip_preprocessing
                disp('Binary image found. Skipping preprocessing.')
                seg_dna = imdna;
            else
                procimg = ml_preprocess(double(imdna),immask,'ml','yesbgsub', options.threshmeth_dna);
                seg_dna = imfill(procimg>0,'hole');
            end
    end
end

if ~isempty(imcell)
    disp('Preprocessing cell membrane image')
    if isbinary(imcell) && options.skip_preprocessing
        disp('Binary image found. Skipping preprocessing.')
        seg_cell = imcell;
        immask = seg_cell;
    else       
        if length(unique(imcell(:))) > 2
            imcell = imcell.*immask;
            celledge = ml_imedge(imcell,[],'ce');
            cellbodyimg = imfill(celledge,'hole');
        else
            cellbodyimg = imcell;
        end
        
        disp('Finding main body in cell membrane image');
        seg_cell = ml_findmainobj(cellbodyimg);
        
        if ~strcmpi(class(seg_cell),'uint8')
            seg_cell = uint8(seg_cell);
            seg_cell = imadjust(seg_cell,stretchlim(seg_cell),[]);
        end
        
        disp('Updating the mask image');
        immask = double(seg_cell).*double(immask>0);
        immask = gather(immask);
        seg_cell = gather(seg_cell);
    end
else
    seg_cell = [];
end

seg_dna = ml_findmainobj(double(seg_dna).*double(immask>0));
seg_dna = uint8(seg_dna);

is_valid = @(x)(strcmpi(class(x),'double') && max(max(x)) == 1 );

if ~is_valid(seg_dna)
    disp('Image must be cast to double before continuing parameterization.');
    seg_dna = double(seg_dna);
    seg_dna(find(seg_dna~=0)) = 1;
end

if ~is_valid(seg_cell)
    disp('Image must be cast to double before continuing parameterization.');
    seg_cell = double(seg_cell);
    seg_cell(find(seg_cell~=0)) = 1;
end
end