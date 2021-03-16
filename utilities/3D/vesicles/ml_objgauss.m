function [mu,sigma] = ml_objgauss(obj,option)
%TZ_OBJGAUSS Estiamte the Gaussian parameters of an object.
%   MU = ML_OBJGAUSS(OBJ) returns a 1x2 vector which is the mean of the
%   gaussian model of the [object] OBJ.
%   
%   [MU,SIGMA] = ML_OBJGAUSS(OBJ) returns both the mean and the convariance
%   matrix of the model.
%   
%   [MU,SIGMA] = ML_OBJGAUSS(OBJ,OPTION) calculates the covariance by the 
%   specified option OPTION:
%       'full' - full covariace matrix
%       'diag' - diagonal covariance matrix
%       'spherical' - spherical covariance matrix
%
%   See also TZ_GAUSSOBJ

%   18-Jan-2006 Initial write T. Zhao
%   16-Dec-2009 Edit T. Peng: 3D compatibility 
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

if nargin < 1
    error('Exactly 1 argument is required')
end

if ~exist('option','var')
    option = 'full';
end

for d = 1:size(obj,2)-1
    mu(d) = ml_wmoment(obj(:,d),obj(:,end),1);
end

if nargout==2
    if size(obj,2)==3
        switch option
            case 'full'
                sigma(1,1) = ml_wmoment(obj(:,1),obj(:,3),2);
                sigma(2,2) = ml_wmoment(obj(:,2),obj(:,3),2);
                sigma(1,2) = tz_wcorr(obj(:,1),obj(:,2),obj(:,3));
                sigma(2,1) = sigma(1,2);
            case 'diag'
                sigma(1,1) = ml_wmoment(obj(:,1),obj(:,3),2);
                sigma(2,2) = ml_wmoment(obj(:,2),obj(:,3),2);
                sigma(1,2) = 0;
                sigma(2,1) = 0;
            case 'spherical'
                sigma(1,1) = (ml_wmoment(obj(:,1),obj(:,3),2)+ ...
                    ml_wmoment(obj(:,2),obj(:,3),2))/2;
                sigma(2,2) = sigma(1,1);
                sigma(1,2) = 0;
                sigma(2,1) = 0;
        end
    else
        switch option
            case 'full'
                sigma(1,1) = ml_wmoment(obj(:,1),obj(:,4),2);
                sigma(2,2) = ml_wmoment(obj(:,2),obj(:,4),2);
                sigma(3,3) = ml_wmoment(obj(:,3),obj(:,4),2);
                sigma(1,2) = tz_wcorr(obj(:,1),obj(:,2),obj(:,4));
                sigma(2,1) = sigma(1,2);
                sigma(1,3) = tz_wcorr(obj(:,1),obj(:,3),obj(:,4));
                sigma(3,1) = sigma(1,3);
                sigma(2,3) = tz_wcorr(obj(:,2),obj(:,3),obj(:,4));
                sigma(3,2) = sigma(2,3);
            case 'diag'
                sigma(1,1) = ml_wmoment(obj(:,1),obj(:,4),2);
                sigma(2,2) = ml_wmoment(obj(:,2),obj(:,4),2);
                sigma(3,3) = ml_wmoment(obj(:,3),obj(:,4),2);
                sigma(1,2) = 0;
                sigma(2,1) = 0;
                sigma(1,3) = 0;
                sigma(3,1) = 0;
                sigma(2,3) = 0;
                sigma(3,2) = 0;
            case 'spherical'
                sigma(1,1) = (ml_wmoment(obj(:,1),obj(:,4),2)+ ...
                    ml_wmoment(obj(:,2),obj(:,4),2)+ ...
                    ml_wmoment(obj(:,3),obj(:,4),2))/3;
                sigma(2,2) = sigma(1,1);
                sigma(3,3) = sigma(1,1);
                sigma(1,2) = 0;
                sigma(2,1) = 0;
                sigma(1,3) = 0;
                sigma(3,1) = 0;
                sigma(2,3) = 0;
                sigma(3,2) = 0;
        end
    end
end

