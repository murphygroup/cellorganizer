function answer = img2tif( img, filename, compression, index )
% IMG2TIF Saves an image to a multichannel tiff. Input argument
% 'compression' can be 'none', 'lzw' or 'packbits'. Parameter index
% is a boolean flag that allows the user to save the image as an
% indexed tiff.
%
% Examples
% >> img2tif( a, 'a.tif' );
% >> img2tif( b, 'b.tif', 'none' );
% >> img2tif( c, 'c.tif', 'lzw' );

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2008-2018 Murphy Lab
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

% March 7, 2012 R.F. Murphy Change stretching from 2% to min-max
% June 5, 2012 M. Mackie Added parameter to deal with indexed images
% August 4, 2012 D. Sullivan Added logical image "contrast stretching"
% August 6, 2012 I. Cao-Berg Fixed statement where logical images were
%                            stretched to the incorrect value
% January 21, 2013 D. Sullivan no stretching seems to be present, added
%                              min-max stretching
% April 29, 2013 D. Sullivan stretchlim function broken, just scaling max
% to 1 now.

%mmackie 6/5/12
answer=false;
if( nargin < 2)
    disp('img2tif:input not enough, exiting method');
    return
end

if( nargin < 3)
    compression = 'none';
    index = false;
end

if( nargin < 4)
    index = false;
end

%Yajing Tang 7/3/13 Check the inputs' types
if (~isnumeric(img) && ~islogical(img))
    disp('img2tif:not a valid image, exiting method');
    return;
end

if (~isstr(filename) || isempty( filename ) )
    disp('img2tif:filename has to be a non-empty string, exiting method');
    return;
end
if isequal(strfind(filename,'.tif'),[]) || (strfind(filename,'.tif')~=(length(filename)-3 ))
    filename=[filename '.tif'];
end

%D. Sullivan 4/29/13 Don't need to do this anymore.
%devins 8/4/12 added logical contrast stretching
if islogical(img)
    %icaoberg 8/6/2012
    %img = img.*256;
    img = img.*255;
end

%D. Sullivan 1/21/13
%no contrast stretching seems to be present, adding min-max stretching
%D. Sullivan 4/29/13 Stretchlim does not seem to be working properly,
%removing it
% low_high = stretchlim(img(:));
% img = img./low_high(2);
% maxval = double(max(img(:)));
% img = double(img)./maxval;

%it will save an image as a mat file just to check contents
%save( [ filename(1:end-4) '.mat' ], 'img' );

%mmackie 6/5/2012
if index
    % Map values to indices
    img_unique = unique(img(:));
    n_img_unique = length(img_unique);
    if n_img_unique <= 2^8
        img2 = zeros(size(img), 'uint8');
    elseif n_img_unique <= 2^16
        img2 = zeros(size(img), 'uint16');
    else
        warning('Unknown image type, renormalizing and saving as double');
        img2 = img;
        for i = 1:size(img2, 3)
            img2(:,:,i) = double(img2(:,:,i)) / max(reshape(img(:,:,i), 1, []));
        end
    end
    if isinteger(img2)
        for i = 1:n_img_unique
            img2(img == img_unique(i)) = i-1;
        end
    end
    img_unique_colors = hsv2rgb([zeros(n_img_unique, 2), img_unique / max(img_unique)]);
    img = img2;
    
    %save first channel
    imwrite( img(:,:,1), gray(64), filename, 'Resolution', 300, 'Compression', compression);
    %append other channels to the first channel
    if size(img,3) > 1
        for i=2:1:size(img,3)
            imwrite( img(:,:,i), gray(64), filename, 'Resolution', 300, 'WriteMode', 'append', 'Compression', compression );
        end
    end
elseif ndims(img)==4
    %save first channel
    imwrite( squeeze(img(:,:,1,:)), filename, 'Resolution', 300, 'Compression', compression );
    
    %append other channels to the first channel
    if size(img,3) > 1
        for i=2:1:size(img,3)
            imwrite( squeeze(img(:,:,i,:)), filename, 'Resolution', 300, 'WriteMode', 'append', 'Compression', compression );
        end
    end
else
    %save first channel
    %D. Sullivan 5/7/13 check the file type for each channel and convert it
    %where appropriate so that outputs are displayed properly.
    if max(img(:))<=1
        imgsave = double(img(:,:,1));
    elseif max(img(:))<=255
        imgsave=uint8(img(:,:,1));
    elseif max(img(:))<=65536
        imgsave=uint16(img(:,:,1));
    else
        warning('Unknown image type, renormalizing and saving as double')
        imgsave = img(:,:,1)./max(max(img(:,:,1)));
    end
    
    imwrite( imgsave, filename, 'Resolution', 300, 'Compression', compression );
    
    %append other channels to the first channel
    if size(img,3) > 1
        for i=2:1:size(img,3)
            %D. Sullivan 5/7/13 check the file type for each channel and convert it
            %where appropriate so that outputs are displayed properly.
            if max(img(:))<=1
                imgsave = double(img(:,:,i));
            elseif max(img(:))<=255
                imgsave=uint8(img(:,:,i));
            elseif max(img(:))<=65536
                imgsave=uint16(img(:,:,i));
            else
                warning('Unknown image type, renormalizing and saving as double')
                imgsave = img(:,:,i)./max(max(img(:,:,i)));
            end
            
            imwrite( imgsave, filename, 'Resolution', 300, 'WriteMode', 'append', 'Compression', compression );
        end
    end
end

answer=true;
return;
end%img2tif
