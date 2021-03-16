function m = tz_estprojmat(p1,p2)
%TZ_ESTPROJMAT Estimate homogeneous projection matrix.
%   M = TZ_ESTPROJMAT(P1,P2) returns the project matrix for transforming
%   from P1 to P2.

%   12-Sep-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

l = tz_homocoeff(p1,p2);
x = tz_homols(l'*l);
m = reshape(x,size(p1,1)+1,length(x)/(size(p1,1)+1))';
