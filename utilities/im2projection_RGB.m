function img = im2projection_RGB( img_orig, opts )
% IM2GPROJECTION creates a sum or mean projection of the input image
%
% Input:
% img = a 3D binary or realvalued image.
% param = struct with a 'method' field that can be set
% to 'mean' if a mean value projection is desired
% Output:
% out_img = a 2D image that contains a projection
% in each dimension of the original image

% Author: Gregory Johnson
%
% Copyright (C) 2014 Murphy Lab
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

if ~exist('opts', 'var')
    opts = struct();
end

opts = ml_initparam(opts, struct(...
                'method', 'mean', ...
                'justz', false, ...
                'cm', @hsv, ...
                'bg', [0,0,0], ...
                'scale_inten', 1, ...
                'colors', [], ...
                'normalize', true ...
                ));
            

out_img = [];
if isempty( img_orig )
    warning( 'Input argument image is empty' );
    return
end

if iscell(img_orig)
    uimg = 1:length(img_orig);
else
    uimg = unique(img_orig);
    uimg(uimg == 0) = [];
end

if isempty(opts.colors)
    colors = opts.cm(length(uimg));
else
    colors = opts.colors;
end

if iscell(img_orig)    
    [m,n,q] = size(img_orig{1});
else
    [m,n,q] = size(img_orig);
end

imgs{1} = zeros((m+q),(n+q),3);

for i = 1:length(uimg)
    try
        
        if iscell(img_orig)
            img = img_orig{i};
        else
            img = img_orig == uimg(i);
        end
        
        if opts.justz
            q = 0;
        end
        
        out_img = zeros((m+q),(n+q));
        
        imx = [];
        imy = [];
        
        switch opts.method
            case 'sum'
                imz = sum(img,3);

                if ~opts.justz
                    imx = squeeze(sum(img,1));
                    imy = squeeze(sum(img,2));
                end

            case 'mean'
                imz = mean(img,3);

                if ~opts.justz
                    imx = squeeze(mean(img,1));
                    imy = squeeze(mean(img,2));
                end
            case 'max' 
                imz = max(img,[],3);
                
                if ~opts.justz
                    imx = squeeze(max(img,[],1));
                    imy = squeeze(max(img,[],2));
                end
            otherwise
                out_img = [];
                warning('Unknown method');
        end
            
        imz = imz * opts.scale_inten;
        imx = imx * opts.scale_inten;
        imy = imy * opts.scale_inten;
        
        out_img(1:m,q+1:end) = imz;
        out_img(1:m,1:q) = imy;
        out_img(m+1:end,q+1:end) = flipud(imx');
    catch err
        out_img = [];
        warning('Unable to make projections');
        getReport( err )
        return
    end
    outimg = repmat(out_img, [1,1,3]);
    
    outimg(:,:,1) = outimg(:,:,1).*colors(i,1);
    outimg(:,:,2) = outimg(:,:,2).*colors(i,2);
    outimg(:,:,3) = outimg(:,:,3).*colors(i,3);
%     outimg = outimg * 255;
    
    imgs{i} = outimg;
end

imgs = cat(4, imgs{:});
imgs = sum(imgs,4);

%color the background according to param.bg
bg = all(imgs==0,3);

bg_color = repmat(bg, [1,1,3]);
for i = 1:3
    bg_color(:,:,i) = bg_color(:,:,i).*opts.bg(i);
end
imgs(bg_color>0) = bg_color(bg_color>0);

if opts.normalize
    for i = 1:3
        chan = imgs(:,:,i);
        imgs(:,:,i) = chan/max(chan(:));
    end
end

imgs(isnan(imgs)) = 0;


img = imgs;


% img = uint8(imgs./(max(imgs(:)))*255);

