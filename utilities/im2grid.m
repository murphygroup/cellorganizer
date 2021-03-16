function img = im2grid(img, panels, hasLabels)
% IM2GRID reshapes a 3D (grayscale) image into a 2D rectangular grid of the
% z-slices
% 
% Usage:
% img = im2grid(img)
%   Visualizes 3D image as a 2D grid using default panel arrangement and
%   labels slice numbers by default.
% img = im2grid(img, panels)
%   Define custom panel arrangement. PANELS=[M,N] which specifies that the
%   output grid is MxN (by default, PANEL=[M,N] for the minimum M+N which
%   gives M*N >= z-size)
% img = im2grid(img, panels, hasLabels)
%   Set HASLABELS=TRUE to number each panel, HASLABELS=FALSE to hide
%   numbering (default HASLABELS=TRUE).
% 
% Author: Kelvin M. Liu-Huang
% Created: April 26, 2016

img = mat2gray(img);
imsize = size(img);
if imsize(3) > 1
    %Compute the number of panels for best viewing
    root = sqrt(imsize(3));
    if ~isdefined('panels')
        %Default panel arrangement is MxN, where M=round(sqrt(size(img,3)))
        %and N=ceil(sqrt(size(img,3)))
        panels = [round(root), ceil(root)];
    end
    
    %Render slice number in top-left corner (optional)
    if ~isdefined('hasLabels') || hasLabels
        minimsize = min(imsize(1:2));
        for iZ = 1:imsize(3)
            %Render number with font size 1/10th image size, white text,
            %transparent box (box color set to black in case)
            img(:,:,iZ) = rgb2gray(insertText( img(:,:,iZ), ...
                [minimsize/100,minimsize/100], iZ, 'FontSize', ...
                min(round(minimsize/10),72), 'TextColor', [1,1,1], ...
                'Boxcolor', [0,0,0], 'BoxOpacity', 0 ));
        end
    end
    
    %Reshape matrix into desired arrangement
    img = padarray( img, [0, 0, panels(1)*panels(2)-imsize(3)], ...
        cast(0,'like',img), 'post' ); %pad panels so there are MxN of them
    img = reshape( img, [imsize(1), panels(2)*imsize(2), panels(1)] );
    img = permute( img, [1,3,2] );
    img = reshape( img, [panels(1)*imsize(1), panels(2)*imsize(2)] );
else
    warning( 'Image array must 3D' );
    img = [];
end

return
