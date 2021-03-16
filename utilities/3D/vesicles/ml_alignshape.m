function [a,t]=ml_alignshape(sh2,sh1,w)
%ML_ALIGNSHAPE Estimate weighted transformation between shapes.
%   A = ML_ALIGNSHAPE(SH2,SH1) returns the rotation and scaling parameters
%   of transformation from shape SH2 to SH1. SH2 and SH1 should have the 
%   same size.
%   
%   A = ML_ALIGNSHAPE(SH2,SH1,W) takes account of weights for each point.
%
%   [A,T] = ML_ALIGNSHAPE(...) also returns the translation vector.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

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
    error('2 or 3 arguments are required')
end

if nargin<3
    w=[];
end

if isempty(w)
    w=ones(size(sh1,1),1);
end

X1=sum(w.*sh1(:,1));
X2=sum(w.*sh2(:,1));
Y1=sum(w.*sh1(:,2));
Y2=sum(w.*sh2(:,2));
W=sum(w);
Z=sum(w.*(sh2(:,1).^2+sh2(:,2).^2));
C1=sum(w.*sum(sh1.*sh2,2));
C2=sum(w.*(sh1(:,2).*sh2(:,1)-sh1(:,1).*sh2(:,2)));

coefmat=[X2,-Y2,W,0;Y2,X2,0,W;Z,0,X2,Y2;0,Z,-Y2,X2];
sol=inv(coefmat)*[X1,Y1,C1,C2]';
a=sol(1:2)';
t=sol(3:4)';
