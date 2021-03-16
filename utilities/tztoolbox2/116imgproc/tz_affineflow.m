function [coef,scale] = tz_affineflow(img1,img2)
%TZ_AFFINEFLOW Affine flow from two images.
%   COEF = TZ_AFFINEFLOW(IMG1,IMG2)
%   
%   See also

%   24-Oct-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

% [x,y] = meshgrid(1:size(img1,1),1:size(img1,2));
% x=x';
% y=y';
% scale(1) = max(x(:));
% scale(2) = max(y(:));
% x = x/scale(1);
% y = y/scale(2);
% 
% 
% Ix = tz_imderiv(img1,'x')*scale(1);
% Iy = tz_imderiv(img1,'y')*scale(2);
% It = tz_imsub(img2,img1);
% 
% 
% 
% x2 = x.^2;
% y2 = y.^2;
% xy = x.*y;
% Ix2 = Ix.^2;
% Iy2 = Iy.^2;
% IxIy = Ix.*Iy;
% IxIt = Ix.*It;
% IyIt = Iy.*It;
% 
% A = [sum(x2(:).*Ix2(:)) sum(xy(:).*Ix2(:)) sum(x(:).*Ix2(:)) ...
%     sum(x2(:).*IxIy(:)) sum(xy(:).*IxIy(:)) sum(x(:).*IxIy(:)); ...
%     0 sum(y2(:).*Ix2(:)) sum(y(:).*Ix2(:)) sum(xy(:).*IxIy(:)) ...
%     sum(y2(:).*IxIy(:)) sum(y(:).*IxIy(:)); ...
%     0 0 sum(Ix2(:)) sum(x(:).*IxIy(:)) sum(y(:).*IxIy(:)) sum(IxIy(:)); ...
%     0 0 0 sum(x2(:).*Iy2(:)) sum(xy(:).*Iy2(:)) sum(x(:).*Iy2(:)); ...
%     0 0 0 0 sum(y2(:).*Iy2(:)) sum(y(:).*Iy2(:)); ...
%     0 0 0 0 0 sum(Iy2(:))];
% A = A+A'-diag(diag(A));
% 
% b = [sum(x(:).*IxIt(:)) sum(y(:).*IxIt(:)) sum(IxIt(:)) ...
%     sum(x(:).*IyIt(:)) sum(y(:).*IyIt(:)) sum(IyIt(:))]';
% 
% coef = -inv(A)*b;

It = tz_imsub(img2,img1);
It(1,:) = [];
It(:,1) = [];
img1(1,:) = [];
img1(:,1) = [];

scale(1) = size(img1,1);
scale(2) = size(img1,2);
imcoords = tz_imcoords(size(img1),scale);

% Ix = tz_imderiv(img1,'x')*scale(1);
% Iy = tz_imderiv(img1,'y')*scale(2);
[Iy,Ix] = gradient(double(img1),1/scale(2),1/scale(1));

A2 = [imcoords(1,:).*Ix(:)'; imcoords(2,:).*Ix(:)'; Ix(:)'; ...
    imcoords(1,:).*Iy(:)'; imcoords(2,:).*Iy(:)'; Iy(:)']';
coef = -inv(A2'*A2)*A2'*It(:);
