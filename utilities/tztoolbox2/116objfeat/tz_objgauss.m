function [mu,sigma] = tz_objgauss(obj,option)
%TZ_OBJGAUSS Obsolete. See ML_OBJGAUSS.
%   MU = TZ_OBJGAUSS(OBJ) returns a 1x2 vector which is the mean of the
%   gaussian model of the [object] OBJ.
%   
%   [MU,SIGMA] = TZ_OBJGAUSS(OBJ) returns both the mean and the convariance
%   matrix of the model.
%   
%   [MU,SIGMA] = TZ_OBJGAUSS(OBJ,OPTION) calculates the covariance by the 
%   specified option OPTION:
%       'full' - full covariace matrix
%       'diag' - diagonal covariance matrix
%       'spherical' - spherical covariance matrix
%
%   See also TZ_GAUSSOBJ

%   18-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_objgauss','ml_objgauss'));

if nargin < 1
    error('Exactly 1 argument is required')
end

if ~exist('option','var')
    option = 'full';
end

mu(1) = ml_wmoment(obj(:,1),obj(:,3),1);
mu(2) = ml_wmoment(obj(:,2),obj(:,3),1);

if nargout==2 
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
end

