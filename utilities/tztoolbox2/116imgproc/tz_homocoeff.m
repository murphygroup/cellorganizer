function l = tz_homocoeff(p1,p2)
%TZ_HOMOCOEFF Build coefficient matrix for projection matrix.
%   L = TZ_HOMOCOEFF(P1,P2) returns a matrix containing the coefficients
%   for solving the projection matrix.  P1 and P2 should have the same
%   number of columns. It is supposed to be from P1 to P2.

%   12-Sep-2005 Initial write T. Zhao

if nargin < 2
    error('Exactly 2 arguments are required')
end

if(size(p1,2) ~= size(p2,2))
    error('The 2 arguments musth have the same number of columns');
end

nPoints = size(p1,2);
dimPoints = size(p1,1);

homoP1 = tz_homoext(p1);
tempL = reshape([homoP1;homoP1],dimPoints+1,nPoints*2)';

l1 = tempL;
l1(2*(1:nPoints),:) = 0;
l2 = [ l1(end,:); l1(1:end-1,:) ];
l3 = diag(p2(:))*tempL;

l = [l1,l2,-l3];

