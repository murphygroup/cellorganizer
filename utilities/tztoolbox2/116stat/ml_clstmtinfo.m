function [info,confmat] = ml_clstmtinfo(label1,label2)
%ML_CLSTMTINFO Calculate mutual information betwee two clustering results.
%   INFO = ML_CLSTMTINFO(LABEL1,LABEL2) returns the mutual information of
%   the two clustering results, which are represented by the [label
%   vectors] LABEL1 and LABEL2.
%
%   [INFO,CONFMAT] = ML_CLSTMTINFO(LABEL1,LABEL2) also returns mutual confusion
%   matrix.
%   
%   See also

%   12-Sep-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

tb = ml_summarytestclassif(label1,label2,max(label2));
confmat = tb;
tb = tb/sum(tb(:));

px = sum(tb,1);
py = sum(tb,2);
pxpy = py*px;
pxpy(pxpy==0) = 1;
pxy = tb;
pxy(pxy==0) = 1;

tb = tb.*log(pxy./pxpy);
info = sum(tb(:));
