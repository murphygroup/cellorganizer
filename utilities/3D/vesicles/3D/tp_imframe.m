function img = tp_imframe(img,w)
% Add blank frame for 3D image, w is the width of frame

if nargin < 2
    w = [1 1 1];
end

[R,C,H] = size(img);
img = cat(3,zeros(R,C,w(3)),img,zeros(R,C,w(3)));
img = cat(1,zeros(w(1),C,H+2*w(3)),img,zeros(w(1),C,H+2*w(3)));
img = cat(2,zeros(R+2*w(1),w(2),H+2*w(3)),img,...
    zeros(R+2*w(1),w(2),H+2*w(3)));