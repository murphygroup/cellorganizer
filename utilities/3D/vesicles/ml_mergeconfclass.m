function confmat = ml_mergeconfclass(confmat,param)
%ML_MERGECONFCLASS Merge classes for a confusion matrix.
%   CONFMAT2 = ML_MERGECONFCLASS(CONFMAT,PARAM) returns a new confusion
%   matrix that is modified from the input confusion matrix CONFMAT by
%   merging classes. PARAM is a structure to specify how to merge the
%   classes. It has the following fields:
%       'colidx' - class indices to merge columns. It could be a row vector
%           or cell array of row vectors. All indices should be different.
%       'rowidx' - class indices to merge rows. It must have the same data
%           lengt with 'colidx'. If it does not exist, it will be set to
%           'colidx'
%       'balance' - 1 for balance the numbers after merging (default), 0
%           not. 
%   
%   See also

%   18-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
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

param = ml_initparam(param,struct('balance',1,'rowidx',{param.colidx}));

if isnumeric(param.colidx)
    param.colidx = {param.colidx};
end

if isnumeric(param.rowidx)
    param.rowidx = {param.rowidx};
end

rmrowidx = [];
rmcolidx = [];

%merge columns
for i=1:length(param.rowidx)
    colidx = param.colidx{i};
    if length(colidx)>1
        confmat(:,colidx(1)) = sum(confmat(:,colidx),2);
        rmcolidx = [rmcolidx colidx(2:end)];
    end
end

%merge rows
for i=1:length(param.rowidx)
    rowidx = param.rowidx{i};
    if length(rowidx)>1
        confmat(rowidx(1),:) = sum(confmat(rowidx,:),1);
        rmrowidx = [rmrowidx rowidx(2:end)];
        if param.balance==1
            confmat(rowidx(1),:) = confmat(rowidx(1),:)/length(rowidx);
        end
    end
end

%remove redundant rows and columns
confmat(:,rmcolidx) = [];
confmat(rmrowidx,:) = [];


