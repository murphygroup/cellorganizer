function seg = region_seg( I, init_mask, max_its, alpha, display, quit_tol, ...
    param )
% Region Based Active Contour Segmentation
%
% seg = region_seg(I,init_mask,max_its,alpha,display,quit_tol)
%
% Inputs: I           2D image
%         init_mask   Initialization (1 = foreground, 0 = bg)
%         max_its     Number of iterations to run segmentation for
%         alpha       (optional) Weight of smoothing term
%                       higer = smoother.  default = 0.2
%         display     (optional) displays intermediate outputs
%                       default = true
%         quit_tol    (optional) quit if fraction of all points that
%                       move between 10 iterations is <= this
%         dirFlag     Direction of propagation. False - outside-in, True -
%                       inside-out
%
% Outputs: seg        Final segmentation mask (1=fg, 0=bg)
%
% Description: This code implements the paper: "Active Contours Without
% Edges" By Chan Vese. This is a nice way to segment images whose
% foregrounds and backgrounds are statistically different and homogeneous.
%
% Example:
% img = imread('tire.tif');
% m = zeros(size(img));
% m(33:33+117,44:44+128) = 1;
% seg = region_seg(img,m,500);

% Author: Shawn Lankton (www.shawnlankton.com)
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

% March 15, 2012 R.F. Murphy Added convergence tolerance
%
% July 1, 2012 G.J. Johnson, Allows for 3D images
%
% March 15, 2012 R.F. Murphy Added convergence tolerance
%
% April 24, 2014 I. Cao-Berg Created PLOT_SWITCH that turns on/off display
% of results. Default is off and this is a hidden parameter only meant for
% developers use.
%
% March 11, 2017 R.F. Murphy Fix display logic

if nargin == 6
    param = [];
    param = ml_initparam( param, struct( 'display', false ) );
    param = ml_initparam( param, struct( 'debug', false ) );
    param = ml_initparam( param, struct( 'verbose', false ) );
end

%-- default value for parameter alpha is .1
if(~exist('alpha','var'))
    alpha = .2;
end

%-- display intermediate outputs if flag set by
%   either argument or param structure
if(~exist('display','var'))
    display = param.display;
end
display = display | param.display;
% PLOT_SWITCH = false;

%-- default behavior is to quit if no points moved between 10 iterations
if(~exist('quit_tol','var'))
    quit_tol = 0;
end
%-- ensures image is 2D double matrix
I = im2graydouble(I);

%-- Create a signed distance map (SDF) from mask
phi = mask2phi(init_mask);
oldphi = phi; % copy for testing exit condition

%-- show starting mask
%icaoberg when these two are set it will save the display to disk. much
%faster than making an actual display
if display
    try
        showCurveAndPhi( I, phi, 0 );
    catch
        disp( 'Unable to open display. Skipping saving of intermediate results plot.' );
    end
end

%icaoberg 4/24/2014
ITERATION_DISPLAY_NUMBER = 100;
%--main loop
if param.verbose
    disp( ['Maximum number of iterations set to: ' num2str(max_its)] );
end
times = 0;
for its = 1:max_its
    tic
    if(mod(its,ITERATION_DISPLAY_NUMBER) == 0)
        if param.verbose
            disp( ['Iteration: ' num2str(its) ] );
        end
    end
    
    idx = find(phi <= 1.2 & phi >= -1.2);  %get the curve's narrow band
    
    if isempty(idx)
         break;
    end
    %-- find interior and exterior mean
    upts = find(phi<=0);                 % interior points
    vpts = find(phi>0);                  % exterior points
    u = sum(I(upts))/(length(upts)+eps); % interior mean
    v = sum(I(vpts))/(length(vpts)+eps); % exterior mean
    
    F = (I(idx)-u).^2-(I(idx)-v).^2;         % force from image information
    curvature = get_curvature(phi,idx);  % force from curvature penalty
    
    dphidt = F./max(abs(F)) + alpha*curvature;  % gradient descent to minimize energy
    
    %-- maintain the CFL condition
    dt = .45/(max(dphidt)+eps);
    
    %-- evolve the curve
    phi(idx) = phi(idx) + dt.*dphidt;
    
    %-- Keep SDF smooth
    phi = sussman(phi, .5);
    
    %-- intermediate output
    if display && (mod(its,ITERATION_DISPLAY_NUMBER) == 0)
        try
            figure(its);
            showCurveAndPhi(I,phi,its);
        catch
            disp( 'Unable to open display. Skipping saving of intermediate results plot.' );
        end
    end
    
    if(~mod(its,10))
        pointsmoved = sum(sum(xor(oldphi<=0,phi<=0)));
        if pointsmoved <= quit_tol*prod(size(I)) break, end % test for convergence
        oldphi = phi; % copy for testing exit condition
    end
    
    times = times + toc;
    if(mod(its,ITERATION_DISPLAY_NUMBER) == 0)
        times = times + toc;
        if param.verbose
            disp( ['Elapsed time: ' num2str(times)] );
        end
        times = 0;
    end
    finalits=its;
end

%-- final output
if display
    try
        showCurveAndPhi(I,phi,its);
        close
    catch err
        getReport( err, 'extended' );
    end
end

%-- make final mask from SDF
seg = phi<=0; %-- Get mask from levelset




%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------


%-- Displays the image with curve superimposed
function showCurveAndPhi(I, phi, i)
mid = round(size(I,3)/2);

% temp_images_directory = [ pwd filesep 'showCurveAndPhi' ];
% if ~exist( temp_images_directory )
%     mkdir( temp_images_directory );
% end

% figure( 'Visible', 'off' );
imshow(I(:,:,mid),[],'initialmagnification',200); hold on;

contour(phi(:,:,mid), [0 0], 'g','LineWidth',4);
contour(phi(:,:,mid), [0 0], 'k','LineWidth',2);
hold off; title([num2str(i) ' Iterations']); drawnow;
% saveas( gcf, [ temp_images_directory filesep ...
%     'iteration' num2str(sprintf('%05d',i)) ], 'png' );

%-- converts a mask to a SDF
function phi = mask2phi(init_a)
phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;

%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
[dimy, dimx] = size(phi);
[y x] = ind2sub([dimy,dimx],idx);  % get subscripts

%-- get subscripts of neighbors
ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

%-- bounds checking
ym1(ym1<1) = 1; xm1(xm1<1) = 1;
yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;

%-- get indexes for 8 neighbors
idup = sub2ind(size(phi),yp1,x);
iddn = sub2ind(size(phi),ym1,x);
idlt = sub2ind(size(phi),y,xm1);
idrt = sub2ind(size(phi),y,xp1);
idul = sub2ind(size(phi),yp1,xm1);
idur = sub2ind(size(phi),yp1,xp1);
iddl = sub2ind(size(phi),ym1,xm1);
iddr = sub2ind(size(phi),ym1,xp1);

%-- get central derivatives of SDF at x,y
phi_x  = -phi(idlt)+phi(idrt);
phi_y  = -phi(iddn)+phi(idup);
phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
    +0.25*phi(iddr)+0.25*phi(idul);
phi_x2 = phi_x.^2;
phi_y2 = phi_y.^2;

%-- compute curvature (Kappa)
curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
    (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);

%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)
[dimy, dimx, c] = size(img);
if(isfloat(img)) % image is a double
    if(c==3)
        img = rgb2gray(uint8(img));
    end
else           % image is a int
    if(c==3)
        img = rgb2gray(img);
    end
    img = double(img);
end

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
% forward/backward differences
a = D - shiftR(D); % backward
b = shiftL(D) - D; % forward
c = D - shiftD(D); % backward
d = shiftU(D) - D; % forward

a_p = a;  a_n = a; % a+ and a-
b_p = b;  b_n = b;
c_p = c;  c_n = c;
d_p = d;  d_n = d;

a_p(a < 0) = 0;
a_n(a > 0) = 0;
b_p(b < 0) = 0;
b_n(b > 0) = 0;
c_p(c < 0) = 0;
c_n(c > 0) = 0;
d_p(d < 0) = 0;
d_n(d > 0) = 0;

dD = zeros(size(D));
D_neg_ind = find(D < 0);
D_pos_ind = find(D > 0);
dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
    + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
    + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;

D = D - dt .* sussman_sign(D) .* dD;

%-- whole matrix derivatives
function shift = shiftD(M)
shift = permute(shiftR(permute(M,[2,1,3])),[2,1,3]);

function shift = shiftL(M)
shift = [ M(:,2:size(M,2),:) M(:,size(M,2),:) ];

function shift = shiftR(M)
shift = [ M(:,1,:) M(:,1:size(M,2)-1,:) ];

function shift = shiftU(M)
shift = permute(shiftL(permute(M,[2,1,3])),[2,1,3]);

function S = sussman_sign(D)
S = D ./ sqrt(D.^2 + 1);
