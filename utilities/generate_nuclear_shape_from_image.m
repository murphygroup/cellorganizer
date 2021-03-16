function [nucimg, nucsurf, outres ] = generate_nuclear_shape_from_image( image, options )
% GENERATE_NUCLEAR_SHAPE_FROM_IMAGE Generates a nuclear shape from the an
% image. 

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2016 Murphy Lab
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('options','var')
    options = [];    
end

options = ml_initparam(options, struct('verbose', false, ...
                            'debug', false, ...
                            'display', false, ...
                            'resolution', [], ...
                            'synthesis', 'framework'));

options.resolution = ml_initparam(options.resolution, struct(...
                                'objects', [1,1,1], ...
                                'cell', [1,1,1], ...
                                'nucleus', [1,1,1]));
                            
[instance,surfmap] = tp_nucimgfeat(image, ...
                    'cylsurf', options );
                            
% downsampling factor used in generating the models
% f = 5;
%D. Sullivan 2/22/13 changing this to be resolution dependent
%note, objects should have 
%D. Sullivan 3/5/13 f is now 1x3 use to account for resolution adjustment
%in each dimension

%icaoberg 8/7/2013
%this change was done to allow users to synthesize only nuclear shape model
%files
%G. Johson 9/7/2013 
%Bug fixes
if strcmpi(options.synthesis,'all')
    f = options.resolution.cell./options.resolution.objects;
    outres = options.resolution.objects;
outres = options.resolution.objects;
elseif strcmpi( options.synthesis, 'framework' )
    f = [1,1,1];
    outres = options.resolution.cell;
else
    f = [1,1,1];
    outres = options.resolution.nucleus;
end

instance.coefs = shiftdim(instance.coefs,-1);

%D. Sullivan 3/5/13 f is now 1x3 use f(3) for z resolution adjustment
%icaober 7/8/2013
if strcmpi( options.synthesis, 'framework' ) || ...
        strcmpi( options.synthesis, 'all' )
    H = ceil(f(3)*instance.height);
else
    H = instance.height;
end

%icaoberg 8/2/2012
%param = ml_initparam(param,...
%    struct('xsize',1024,'ysize',1024,'zsize',H+1,'samp_rate',360,'debug',true));
factor = 1.25;
% factor = 1;
%D. Sullivan 7/6/13 - this should not be hardcoded at factor*1024. it makes
%no sense. 
options = ml_initparam(options,...
    struct('xsize',factor*1024,'ysize',factor*1024,'zsize',H+1,'samp_rate',360,'debug',true));

xcenter = options.xsize / 2;
ycenter = options.ysize / 2;

delta = 2*pi/options.samp_rate;
Phi = -pi:delta:pi;

% generate one nucleus surface at the specified height
Z = 0:(1/H):1;
[Phi_grid, Z_grid] = meshgrid(Phi,Z);
mesh_data = [Z_grid(:), Phi_grid(:)]';
nucsurf = f(1) * reshape(fnval(instance,mesh_data),[length(Z),length(Phi)]);

%icaoberg 7/1/2013
if options.debug && options.display
    plotcylsurf(nucsurf,delta);
end

% generate another one at the mean height for later use in finding the
% nuclear position
% MEAN_HEIGHT = 85;
% Z = 0:(1/MEAN_HEIGHT):1;
% [Phi_grid, Z_grid] = meshgrid(Phi,Z);
% mesh_data = [Z_grid(:), Phi_grid(:)]';
%D. Sullivan 3/5/13 f is now 1x3 use f(1) for x and y resolution adjustment
% nucsurf2 = f(1) * reshape(fnval(instance,mesh_data),[length(Z),length(Phi)]);
%nucsurf2 = nucsurf2/2;

nucimg = uint8(zeros(options.ysize,options.xsize,options.zsize));
for i = 1:size(nucsurf,1)
    [x,y] = pol2cart(Phi,nucsurf(i,:));
    x = x + xcenter;
    y = y + ycenter;
    sliceimg = zeros(options.ysize,options.xsize);
    for t = 1:length(x)-1
        rpts = round(linspace(y(t),y(t+1),100));
        cpts = round(linspace(x(t),x(t+1),100));
        index = sub2ind([options.ysize,options.xsize],rpts,cpts);
        sliceimg(index) = 255;
    end
    nucimg(:,:,i) = imfill(sliceimg,'holes');
end

%icaoberg 7/1/2013
if options.debug && options.display
    try
        for i=1:size(nucimg,3)
            imshow(nucimg(:,:,i),[0 255]);
            pause(0.1)
        end
    catch
        disp('Unable to display image');
    end
end

%blankimg = zeros(param.ysize,param.xsize);
%nucimg = cat(3,blankimg,nucimg,blankimg);
