function map = tp_surfmap(img,delta)
% TP_SURFMAP convert a 3D surface into a 2D map using Cartesian-Cylindrical
% coordinate transformation. DELTA is the resampling step size of theta

img = img > 0;
[H,W,S] = size(img);
ctrX = H/2;
ctrY = W/2;

Phi = -pi:delta:pi;
N = length(Phi);
R = zeros(S,N);

for k = 1:S
    imgk = img(:,:,k);
    % Dectect boundary
    bounds = bwboundaries(imgk,8);
    %Serena 10/20
    if isempty(bounds)
        continue
    end
    bound = bounds{1};
    boundX = bound(:,1) - ctrX;
    boundY = bound(:,2) - ctrY;
    % Resampling
    [theta,rho] = cart2pol(boundX,boundY);
    [theta,uidx] = unique(theta);
    R(k,:) = interp1([theta-2*pi;theta;theta+2*pi],...
        repmat(rho(uidx),[3,1]),Phi);
end

map = R;
