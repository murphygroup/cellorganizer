function [b,res]=tz_decomobjcom2(y,x1,x2,norm)
%TZ_DECOMOBJCOM2 Use regression method to solve multinomial mixture.
%   B = TZ_DECOMOBJCOM2(Y,X1,X2,NORM)
%   
%   [B,RES] = TZ_DECOMOBJCOM2(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [b,res]=tz_decomobjcom2(y,x1,x2,norm)
%
%OVERVIEW:
%   use regression method to solve multinomial mixture
%PARAMETERS:
%   y - mixture data
%   x1 - component data 1
%   x2 - component data 2
%   norm - normalize x1 and x2 or not
%RETURN:
%   b - coefficients
%   res - residual
%Created by tingz on Mar.3, 2004

if ~exist('norm','var')
    norm=0;
end

if length(y)==1
    warning('There are only one kind of object.')
    return
end

if norm~=0
    normy=tz_normobjcom(y);
    normx1=tz_normobjcom(x1);
    normx2=tz_normobjcom(x2);
    z=normy-normx2;
    x=normx1-normx2;
else
    z=y-x2;
    x=x1-x2;
end

[b,bint,res]=regress(z',x',0.05);

if b>1
    b=1;
end

if b<0
    b=0;
end
