function [spfeat,surfmap] = tp_nucimgfeat( segdna,...
    method, options )
%TP_NUCIMGFEAT Calculates spline features on a nuclear image

% Tao Peng
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
%`
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

% March 22, 2012 Added param.downsample so that user can select a
%                downsampling scale. The default is [5 5 1] which is
%                the one in the original generative model paper
%
% March 23, 2012 Fixed a bug when getting cellnum from imgdir
%
% March 26, 2012 Added image rotation after preprocessing and well as
%                removal of blank slices on top and bottom of preprocessed
%                image
%
% July 23, 2012 R.F. Murphy Add debug code to show steps and final spline
%
% July 26, 2012 I. Cao-Berg If the number of slices in the nuclear image that have
%                           fluorescence is less than 4, then ignore such image
% Jan 1, 2013 I. Cao-Berg Plots will be made iff display flag is set to true
%
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
%
% June 7-13 2013 D. Sullivan changed to per-cell type computations. Major
%                            changes include passing a presegmented file
%                            path rather than raw image files and adding a
%                            'currfile' parameter to track what loop
%                            iteration we are currently on in parallel.
%                            Also passing in a tempdirectory
%                            structure as part of param.
%
% Aug 12, 2013 G. Johnson    Corrected save file name. Improved for
%                            readability
%
% March 12, 2017 R.F. Murphy Pause and then close display when done

if ~exist('options', 'var') || isempty(options)
    options = [];
end

options = ml_initparam(options, struct('verbose', false, ...
    'debug', false, ...
    'display', false, ...
    'downsampling', [1,1,1]));

%icaoberg - march 26, 2012 added image rotation
img = double(segdna);

%icaoberg 26/7/2012
nnzpixel = squeeze(sum(sum(segdna)));

if length(find(nnzpixel)>0) < 4
    spfeat = struct([]);
    surfmap = [];
    return
end

%icaoberg 06/02/2013
if options.debug && options.display
    try
        subplot(2,2,1);
        for i=1:size(img,3)
            try
                imshow(img(:,:,i),[0 1])
                pause(0.1);
            catch
                disp('Unable to open display.');
            end
        end
    catch
    end
end

try
    % image rotation
    [maxfluo,supid] = max(nnzpixel);
    supshape = img(:,:,supid);
    theta = ml_majorangle(supshape)*180/pi;
    if options.verbose
        disp(['Major angle for XY rotation is ' num2str(theta) '.']);
    end
    
    %icaoberg 5/15/2013
    %had to change this snippet because ml_rotate might change the size if
    %the rotated image
    img2 = [];
    for k = find(nnzpixel,1):find(nnzpixel,1,'last')
        temp = ml_rotate(img(:,:,k),-theta);
        img2 = cat( 3, img2, temp );
    end
    img = img2;
    clear img2;
    
    %icaoberg 1/30/2013
    if options.debug && options.display
        try
            subplot(2,2,2);
            for i=1:size(img,3)
                imshow(img(:,:,i),[0 1])
                pause(0.1);
            end
        catch
            disp('Unable to open display.');
        end
    end
catch error
    spfeat = [];
    if options.verbose
        disp( ['Could not rotate image.']);
    end
    getReport( err )
    return
end

% spline surface fitting
img = tp_imtight(img);
switch method
    case 'medsurf'
        [medplane,height,mask] = tp_imaxisplane(img);
        [medfeat,heifeat] = tp_medplanefeat(medplane,height,mask);
        spfeat = struct('medsurf',medfeat,'height',heifeat);
    case 'cylsurf'
        delta = 2*pi/360.;
        surfmap = tp_surfmap(img,delta);
        if options.debug && options.display
            try
                subplot(2,2,3);
                plotcylsurf(surfmap,delta);
            catch
                disp( 'Unable to open display.');
            end
        end
        spfeat = tp_spcylfeat(surfmap);
        
        %icaoberg 1/30/2013
        if options.debug && options.display
            try
                subplot(2,2,4);
                Phi = -pi:delta:pi;
                H = size(surfmap,1); Z = 0:(1/H):1;
                [Phi_grid, Z_grid] = meshgrid(Phi,Z);
                mesh_data = [Z_grid(:), Phi_grid(:)]';
                nucsurf = reshape(fnval(spfeat,mesh_data),[length(Z),length(Phi)]);
                plotcylsurf(nucsurf,delta);
                pause(10)
                close
            catch
                warning( 'Unable to open display.' );
            end
        end
end


end
