%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 04/15/2018 xruan change the input and output, and not save intermediate
% results. 
% 04/27/2018 xruan change the rotation only along z-axis, which is more
% suitable for cell images, and also rotate in the center of image.
% 04/30/2018 add the check of skewness to make sure it is not flipped.
% 10/09/2018 add the option to provide rotation matrix for object and only align the parameters. 
% 02/24/2019 add options for rotation either in xy-plane (default), xz, yz or xyz 

function [vertices, sph_verts, faces, fvec, fvec_obj, R, centers_1]=align_FOE_2(faces, vertices, sph_verts, fvec, options)

% load(filename);
% [path,name,ext] = fileparts(filename);
% 
% if ~exist('faces', 'var') | ~exist('vertices', 'var') | ~exist('sph_verts', 'var') | ~exist('fvec', 'var')
%     disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
%     return;
% end

switch deblank(char(options.CPoint))
    case 'x'
        blue = 1;
    case 'y'
        blue = 2;        
    case 'z'
        blue = 3;        
end

switch deblank(char(options.NPole))
    case 'x'
        yellow = 1;
    case 'y'
        yellow = 2;        
    case 'z'
        yellow = 3;        
end

blueyellow = [blue yellow];
degree = options.MaxSPHARMDegree;


% dirName = [options.OutDirectory '/alignParam'];
% if ~exist(dirName,'dir')
%     mkdir(dirName);
% end

% new_name = sprintf('%s/%sFOE_prm.mat',dirName,name(1:end-3));
% if exist(new_name,'file')
%     prompt = {'Enter new filename:'};
%     dlg_title = 'New File Name';
%     num_lines = 1;
%     def = {new_name};
%     answer = inputdlg(prompt,dlg_title,num_lines,def);    
%     new_name = answer{1};
% end
% save(new_name, 'vertices', 'sph_verts', 'faces', 'fvec','expts');


% rotate in the object space
% svs = [0 0 1; 1 0 0; 0 0 -1]; % north pole, intersection of dateline and equator, south pole
% R = object_rotate_R(vs);
if isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix && isfield(options, 'R')
	R = options.R;
else 
    % rotate the parameter space
    disp('<< Rotate the parameter space >>');
    [fvec, sph_verts, expts] = param_rotate(fvec, vertices, sph_verts, faces, degree, blueyellow);
    svs = [1 0 0; 0 0 1; -1 0 0]; % north pole, intersection of dateline and equator, south pole
    Z = calculate_SPHARM_basis(svs, 1);

    disp('<< Rotate the object space >>');
    vs = real(Z(:,2:4)*fvec(2:4,:));
    
    switch options.rotation_plane
        case {'xy', 'yz', 'xz'}
            [R] = object_rotate_R_z(vs, options.rotation_plane);
        case 'xyz'
            R = object_rotate_R(vs);
    end
end
% R = object_rotate_R([ellipAxes(:,1)';ellipAxes(:,3)']);

if isfield(options, 'use_given_center') && options.use_given_center
    center_obj = options.center;
else 
	center_obj = real(fvec(1, :)) ./ sqrt(4 * pi);
end 

fvec_obj = fvec*R'; 
fvec_obj(1, :) = center_obj .* sqrt(4 * pi);
centers_1 = center_obj;
vertices_c = vertices - centers_1;
vertices_c = vertices_c*R'; 
vertices = vertices_c + centers_1;

if isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix 
    [fvec, sph_verts, expts] = param_rotate(fvec, vertices, sph_verts, faces, degree, blueyellow);
    [fvec_obj, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');        
end


if false
    [vs, fs]=SpiralSampleSphere(4002);
    deg = 31;
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_obj);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    xlabel('x')    
    ylabel('y')
    title('original')
end


% check whether need to flip the meshes
if ~false && ~(isfield(options, 'use_given_rotation_matrix') && options.use_given_rotation_matrix)
    skew = skewness(vertices);

    flipdim = skew(1) < 0;
    if any(flipdim)
        R_skew = diag([-1, -1, 1]);
        vertices_c = vertices - centers_1;
        vertices_c = vertices_c*R_skew'; 
        vertices = vertices_c + centers_1;  
        R = R_skew * R;
        % sph_verts = sph_verts * diag([1, 1, -1]);
        [fvec, sph_verts, expts] = param_rotate(fvec, vertices, sph_verts, faces, degree, blueyellow);
        [fvec_obj, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');        
    end
    
    % 02/24/2019 xruan also check the skewness of z
    if strcmp(options.rotation_plane, 'xyz') && skew(3) < 0 
        R_skew = diag([1, -1, -1]);
        % fvec_obj = fvec_obj * R_skew';
        vertices_c = vertices - centers_1;
        vertices_c = vertices_c*R_skew'; 
        vertices = vertices_c + centers_1;  
        R = R_skew * R;    
        % sph_verts = sph_verts * diag([1, 1, -1]);        
        [fvec, sph_verts, expts] = param_rotate(fvec, vertices, sph_verts, faces, degree, blueyellow);
        [fvec_obj, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');                
    end    
end

if false
    [vs, fs]=SpiralSampleSphere(4002);
    deg = 31;
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_obj);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    xlabel('x')    
    ylabel('y')    
end

if false
    % sph_verts = sph_verts * diag([1, 1, -1]);
    sph_verts = sph_verts * diag([1, 1, -1]);    
    [fvec_obj_1, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');
    [vs, fs]=SpiralSampleSphere(4002);
    deg = 31;
    Zs = calculate_SPHARM_basis(vs, deg);
    Zvert_pdm = real(Zs*fvec_obj_1);
    figure, patch('vertices', Zvert_pdm, 'faces', fs, 'FaceVertexCData',jet(size(Zvert_pdm,1)),'FaceColor','interp');
    xlabel('x')    
    ylabel('y')    
end


return;

%
% Factoring a Rotation Matrix as Rz*Ry*Rx (counterclockwise when looking towards the origin)
%

function [thetaX, thetaY, thetaZ] = factor_rot_xyz(R)

thetaY = asin(R(1,3));
if (thetaY < pi/2)
    if (thetaY > -pi/2)
        thetaX = atan2(-R(2,3),R(3,3));
        thetaZ = atan2(-R(1,2),R(1,1));
    else
        disp('WARNING (Factor Rotation): thetaY = -pi/2, not a unique solution, set thetaZ = 0');
        thetaX = -atan2(R(2,1),R(2,2));
        thetaZ = 0;
    end
else
    disp('WARNING (Factor Rotation): thetaY = pi/2, not a unique solution, set thetaZ = 0');
    thetaX = atan2(R(2,1),R(2,2));
    thetaZ = 0;
end

% disp(sprintf('Factor rotation xyz: %0.2f %0.2f %0.2f',thetaX/pi,thetaY/pi,thetaZ/pi));

return;

%
% rotation matrix in object space
%

function R = object_rotate_R(vs)

% fix north pole
[PHI,THETA] = cart2sph(vs(1,1),vs(1,2),vs(1,3));
ind = find(PHI<0); PHI(ind) = PHI(ind)+2*pi;
THETA = pi/2-THETA;
alpha = -PHI; beta = -THETA;
R = rotate_mat(0, beta, 0)*rotate_mat(0, 0, alpha);
vs = vs*R';
% fix intersection;
[PHI,THETA] = cart2sph(vs(2,1),vs(2,2),vs(2,3));
gamma = -PHI;
R1 = rotate_mat(0, 0, gamma); R = R1*R;
vs = vs*R1';

return;


% rotate only along z-axis

function [R] = object_rotate_R_z(vs, rotation_plane)

if nargin < 2
    rotation_plane = 'xy';
end

% fix north pole
[PHI,THETA] = cart2sph(vs(1,1),vs(1,2),vs(1,3));
ind = find(PHI<0); PHI(ind) = PHI(ind)+2*pi;
THETA = pi/2-THETA;
alpha = -PHI; beta = -THETA;
R = rotate_mat(0, beta, 0)*rotate_mat(0, 0, alpha);
vs = vs*R';
% fix intersection;
[PHI,THETA] = cart2sph(vs(2,1),vs(2,2),vs(2,3));
gamma = -PHI;
R1 = rotate_mat(0, 0, gamma); R = R1*R;
vs = vs*R1';

% 02/26/2019 add compatible support for old versions (<9.5) 
% only rotate along a given axis
if verLessThan('matlab','9.5')
    eul = rotm2eul(R, 'ZYX');
    switch rotation_plane
        case 'xy'
            R = eul2rotm([eul(1), 0, 0], 'ZYX');
        case 'yz'
            R = eul2rotm([0, 0, eul(3)], 'ZYX');
        case 'xz'
            R = eul2rotm([0, eul(2), 0], 'ZYX');
    end    
else
    eul = rotm2eul(R, 'XYZ');
    switch rotation_plane
        case 'xy'
            R = eul2rotm([0, 0, eul(3)], 'XYZ');
        case 'yz'
            R = eul2rotm([eul(1), 0, 0], 'XYZ');
        case 'xz'
            R = eul2rotm([0, eul(2), 0], 'XYZ');
    end
end

return;


%
% Parameter space rotation using degree 1 ellipsoid
% 

function [fvec, sph_verts, expts] = param_rotate(fvec, vertices, sph_verts, faces, degree, blueyellow)

% calculate matrix A
coeffs = fvec(2:4,:);
A(:,1) = (coeffs(1,:)-coeffs(3,:))';
A(:,2) = -(coeffs(1,:)+coeffs(3,:))'*1i; % there is a typo in the paper, should be - here.
A(:,3) = sqrt(2)*coeffs(2,:)';
A = real(A*sqrt(3)/(2*sqrt(2*pi)));

% SVD to find rotation and scaling matrics
%   [U,S,V] = svd(X) produces a diagonal matrix S of the same dimension as X, with
%   nonnegative diagonal elements in decreasing order, and unitary matrices U and V so
%   that X = U*S*V'.
% need to rotate object space first.
[U,S,V] = svd(A);
% extremum and saddle points
expts = U*S;

% set up rotation matrix
R = V(:,[3 2 1])';
if (det(R)<0)
    disp('WARNING (Parameter Space): rotoinversion!! Change back to pure rotation');
    R(2,:) = R(2,:)*(-1);
end

[thetaX, thetaY, thetaZ] = factor_rot_xyz(R');

% Question: why multiply thetaY by -1?
R = rotate_mat(thetaX, -thetaY, thetaZ);
fprintf('Rotation xyz: %0.2f %0.2f %0.2f\n',thetaX/pi,thetaY/pi,thetaZ/pi);

% the parameter space rotation
sph_verts = (R*sph_verts')';
% create new spharm descriptor (degree 1 is enough)
[fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');

% calculate the blue point
svs = [0 0 1; 1 0 0]; % yellow (north pole), blue (intersection)
Z = calculate_SPHARM_basis(svs, 1);

vs = real(Z(:,2:4)*fvec(2:4,:));
fprintf('blue: (%f, %f, %f)\n',vs);
% blue point and yellow point should be on the positive side of 
% x (blue=1), y (blue=2), or z (blue=3) axis in the object space
% i.e., vs(blue) and vs(yellow) should be >0
blue = blueyellow(1); yellow = blueyellow(2);
if vs(2,blue)<0 || vs(1,yellow)<0
    if vs(2,blue)<0
        Rfix = rotate_mat(0, 0, pi);
    else
        Rfix = eye(3);
    end
    if vs(1,yellow)<0
        Rfix =  rotate_mat(pi, 0, 0)*Rfix;
    end
	% the parameter space rotation
	sph_verts = (Rfix*sph_verts')';
	% create new spharm descriptor (degree 1 is enough)
    [fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, 1, '', '');
end

[fvec, d, Z_tmp1, name_temp] = create_SPHARM_des_LSF(vertices, [], sph_verts, degree, '', '');

return;


