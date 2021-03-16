function [z,zmean,zsdev] = ml_zscore(d,dmean,dsdev)

%ML_ZSCORE normalizes the data into 0 mean and 1 dev
%   [Z,ZMEAN,ZSDEV] = ML_ZSCORE(X) normalizes each column
%   of X into zero meand and unit deviance. The normalized
%   data returns as Z. ZMEAN and ZSDEV are original means
%   and deviances of X.
%   [Z,ZMEAN,ZSDEV] = ML_ZSCORE(X,DMEAN) normalizes the data
%   by specified DMEAN.
%   [Z,ZMEAN,ZSDEV] = ML_ZSCORE(X,DMEAN,DSDEV) normalizes the data
%   by specified DMEAN and DSDEV.
   
% Copyright (C) 2006  Murphy Lab
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

%HISTORY
%   16-May-2005 Initial write TINGZ

n = size(d,1) ;

if ~exist('dmean','var')
    dmean=[];
end
if isempty(dmean);
    zmean=mean(d,1);
else
    zmean=dmean;
end

if ~exist('dsdev','var')
    dsdev=[];
end
if isempty(dsdev)
    zsdev=std(d,0,1);
    %check constant index
    constidx=ml_constidx(d);
    if ~isempty(constidx)
        zsdev(constidx)=1;
    end
else
    zsdev=dsdev;
end

z = (d-(ones(n,1)*zmean)) ./ (ones(n,1)*zsdev);
