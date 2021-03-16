function result = tp_gennucshape( instance,options )
% TP_GENNUCSHAPE generates a nuclear shape from the tensor-product spline
% surface structure

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2012-2019 Murphy Lab
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

% July 23, 2012 R.F. Murphy Fix bug that caused approximately only the
%                   bottom half of nucleus to be returned
% August 1, 2012 R.F.Murphy Don't shrink nuclear size - throws off scaling
%                   of protein objects relative to nucleus and cell
% Feb 22, 2013 D. Sullivan made scaling factor resolution dependent
% March 5, 2013 D. Sullivan f is now 1x3 use to account for resolution
%               adjustment in each dimension


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

nuc_height = length(Z);
% Create mesh
n_angles = length(Phi)-1;
n_slices = nuc_height;
nucmesh = struct('vertices', zeros(n_slices * n_angles + 2, 3), 'faces', zeros((n_slices-1) * n_angles * 2 + n_angles * 2, 3));
% Bottom
bottom_face_indices = (1:n_angles) + (n_slices-1) * n_angles * 2;
bottom_center_vertex_index = size(nucmesh.vertices, 1) - 1;
nucmesh.faces(bottom_face_indices, 1) = 1:n_angles;
nucmesh.faces(bottom_face_indices, 2) = circshift(nucmesh.faces(bottom_face_indices, 1), 1);
nucmesh.faces(bottom_face_indices, 3) = bottom_center_vertex_index;
nucmesh.vertices(bottom_center_vertex_index, 1) = xcenter;
nucmesh.vertices(bottom_center_vertex_index, 2) = ycenter;
nucmesh.vertices(bottom_center_vertex_index, 3) = 1;
% Top
top_face_indices = (1:n_angles) + (n_slices-1) * n_angles * 2 + n_angles;
top_center_vertex_index = size(nucmesh.vertices, 1);
nucmesh.faces(top_face_indices, 1) = (1:n_angles) + (n_slices-1) * n_angles;
nucmesh.faces(top_face_indices, 2) = circshift(nucmesh.faces(top_face_indices, 1), -1);
nucmesh.faces(top_face_indices, 3) = top_center_vertex_index;
nucmesh.vertices(top_center_vertex_index, 1) = xcenter;
nucmesh.vertices(top_center_vertex_index, 2) = ycenter;
nucmesh.vertices(top_center_vertex_index, 3) = n_slices;
% Sides
for i = 1:n_slices-1
    side_face_indices1 = (1:n_angles) + (i-1) * n_angles * 2;
    nucmesh.faces(side_face_indices1, 1) = (1:n_angles) + (i-1) * n_angles;
    nucmesh.faces(side_face_indices1, 2) = circshift(nucmesh.faces(side_face_indices1, 1), -1);
    nucmesh.faces(side_face_indices1, 3) = (1:n_angles) + i * n_angles;
    side_face_indices2 = side_face_indices1 + n_angles;
    nucmesh.faces(side_face_indices2, 1) = nucmesh.faces(side_face_indices1, 1);
    nucmesh.faces(side_face_indices2, 2) = nucmesh.faces(side_face_indices1, 3);
    nucmesh.faces(side_face_indices2, 3) = circshift(nucmesh.faces(side_face_indices2, 2), 1);
end

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
any_points_outside_image = false;
for i = 1:size(nucsurf,1)
    [x,y] = pol2cart(Phi,nucsurf(i,:));
    
    x = x + xcenter;
    y = y + ycenter;
    
    vertices_indices = (1:n_angles) + (i-1) * n_angles;
    nucmesh.vertices(vertices_indices, 1:2) = [x(1:end-1)', y(1:end-1)'];
    nucmesh.vertices(vertices_indices, 3) = i;
    
    sliceimg = zeros(options.ysize,options.xsize);
    for t = 1:length(x)-1
        rpts = round(linspace(y(t),y(t+1),100));
        cpts = round(linspace(x(t),x(t+1),100));

        % T. Buck 2019-01-29 Prevent "Out of range subscript." error but keep filling functionality
        rpts = max(rpts, 1);
        rpts = min(rpts, options.ysize);
        cpts = max(cpts, 1);
        cpts = min(cpts, options.xsize);

        index = sub2ind([options.ysize,options.xsize],rpts,cpts);
        sliceimg(index) = 255;
    end
    nucimg(:,:,i) = imfill(sliceimg,'holes');
end
if any_points_outside_image
    warning('CellOrganizer: Out of range subscripts when adding pixels to image')
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


result = struct('nucimg', nucimg, 'nucsurf', nucsurf, 'outres', outres, 'options', options, 'nucmesh', nucmesh);
