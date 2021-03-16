function pts = ml_mxs2crd(medaxis,width)
%ML_MXS2CRD Convert medial axis shape to coordinates
%   PTS = ML_MXS2CRD(MEDAXIS,WIDTH) returns an array of points which are
%   represented by the medial axis shape with medial axsi MEDAXIS ans width
%   WIDTH.
%   
%   See also

%   30-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

% Copyright (C) 2007  Murphy Lab
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

% pts1(:,1) = medaxis(:,1);
% pts1(:,2) = medaxis(:,2)-floor(width'/2); % **^*
% 
% pts2(:,1) = medaxis(:,1);
% pts2(:,2) = medaxis(:,2)+floor(width'/2-0.5); % **^*
% 
% pts = [pts1;flipud(pts2)];
% pts(end+1,:) = pts(1,:);

pts1(:,1) = medaxis(:,1);
pts1(:,2) = medaxis(:,2)-floor(width'/2); % **^*

pts2(:,1) = medaxis(:,1);
pts2(:,2) = medaxis(:,2)+floor(width'/2-0.5); % **^*

pts12 = [pts1;flipud(pts2)];
x = pts12(:,2); % x-coord
y = pts12(:,1); % y-coord
hp = length(x)/2; % half-point

%%%%%%%%%%%%%%%%
% interp1 method
%%%%%%%%%%%%%%%%

% increment = 0.001;
% spread = 100;
% xit = (x(hp):increment:x(hp+1))'; % input x-coord for top
% xt = x(hp-spread:hp+spread);
% yt = y(hp-spread:hp+spread);
% yit = interp1(xt,yt,xit,'spline'); % outputs y-coord for top
% xib = (x(end):-increment:x(1))'; % input x-coord for bottom
% xb = [x(end-spread:end);x(1:spread)];
% yb = [y(end-spread:end);y(1:spread)];
% yib = interp1(xb,yb,xib,'spline'); % outputs y-coord for bottom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manual interpolation method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Top
% m1T = (y(hp)-y(hp-1))/(x(hp)-x(hp-1)); % top slope 1
% m2T = (y(hp+2)-y(hp+1))/(x(hp+2)-x(hp+1)); % top slope 2
% theta1T = atan(m1T); % top theta 1
% theta2T = atan(m2T); % top theta 2
% dthetaT = 0.0001; % top theta increment
% numpiecesT = length(theta1T:-0.0001:theta2T);
% xdistT = x(hp+1)-x(hp);
% dxT = xdistT/numpiecesT;
% idxT = 0;
% dyT = zeros(1,numpiecesT);
% newxT = zeros(1,numpiecesT);
% for thetaT = theta1T:-0.0001:theta2T
%     idxT = idxT+1;
%     dyT(idxT) = dxT*tan(thetaT);
%     newxT(idxT) = x(hp)+(dxT*idxT);
% end
% newyT = y(hp)+cumsum(dyT);
% 
% % Bottom
% m1B = (y(end)-y(end-1))/(x(end)-x(end-1));
% m2B = (y(2)-y(1))/(x(2)-x(1));
% theta1B = atan(m1B);
% theta2B = atan(m2B);
% dthetaB = dthetaT;
% numpiecesB = length(theta1B:-dthetaB:theta2B);
% xdistB = x(end)-x(1);
% dxB = xdistB/numpiecesB;
% idxB = 0;
% dyB = zeros(1,numpiecesB);
% newxB = zeros(1,numpiecesB);
% for thetaB = theta1B:-dthetaB:theta2B
%     idxB = idxB+1;
%     dyB(idxB) = dxB*tan(thetaB);
%     newxB(idxB) = x(end)-(dxB*idxB);
% end
% newyB = y(end)-cumsum(dyB);
% 
% % Output
% pts = [pts1;[newyT' newxT'];flipud(pts2);[newyB' newxB']];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting to Cubic Function Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xT = x(hp-10:hp+10);
yT = y(hp-10:hp+10);
pT = polyfit(xT,yT,3);
xiT = x(hp):0.001:x(hp+1);
yiT = pT(1).*xiT.^3+pT(2).*xiT.^2+pT(3).*xiT+pT(4);

xB = [x(end-10:end); x(1:10)];
yB = [y(end-10:end); y(1:10)];
pB = polyfit(xB,yB,3);
xiB = x(end):-0.001:x(1);
yiB = pB(1).*xiB.^3+pB(2).*xiB.^2+pB(3).*xiB+pB(4);

pts = [pts1;[yiT' xiT'];flipud(pts2);[yiB' xiB']];


%%%%%%%%%%%%%%%%%%
% Smoothing Method
%%%%%%%%%%%%%%%%%%

% Output
% pts = [pts1;[yit xit];flipud(pts2);[yib xib]];
% xi = pts(:,2);
% yi = pts(:,1);

% window = 200;
% type = 3;
% xst = fastsmooth(xit,window,type,1);
% yst = fastsmooth(yit,window,type,1);
% xsb = fastsmooth(xib,window,type,1);
% ysb = fastsmooth(yib,window,type,1);

% window = 10;
% xst = fastbsmooth(xit,window);
% yst = fastbsmooth(yit,window);
% xsb = fastbsmooth(xib,window);
% ysb = fastbsmooth(yib,window);

% pts = [pts1;[yst xst];flipud(pts2);[ysb xsb]];
% xs = ptsmooth(:,2);
% ys = ptsmooth(:,1);

% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(1,2,1);
% plot(x,y);
% subplot (1,2,2);
% plot(pts(:,2),pts(:,1));
% figure;

