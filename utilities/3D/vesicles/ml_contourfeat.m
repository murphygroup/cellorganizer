function feats = ml_contourfeat(pts)
%ML_CONTOURFEAT Calculate features for a contour.
%   FEATS = ML_CONTOURFEAT(PTS) returns a vector of features of a contour
%   represented by PTS, which is a 2-column matrix. There are 3 features,
%   including 'total length','no. of breaks' and 'no. of edge points'.
%   
%   See also ML_TRACECONTOUR

%   06-Apr-2005 Initial write T. Zhao
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

if nargin < 1
    error('Exactly 1 argument is required')
end

feats(1)=size(pts,1);
pts=[pts;pts(1,:)];

nbdists=max(abs(pts(1:end-1,:)-pts(2:end,:)),[],2);
feats(2)=sum(nbdists>1);

edgelen=[sum(pts(:,1)==0),sum(pts(:,2)==0),...
        sum(pts(:,1)==max(pts(:,1))),sum(pts(:,2)==max(pts(:,2)))];
feats(3)=max(edgelen);