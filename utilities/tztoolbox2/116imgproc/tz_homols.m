function x = tz_homols(v)
%TZ_HOMOLS Homogeneous least squares.
%   X = TZ_HOMOLS(V) returns a vector which is to minimize X'*V*X and
%   has |X|=1.

%   12-Sep-2005 Initial write T. Zhao

if nargin < 1
    error('Exactly 1 argument is required')
end

[eigenVectors,eigenValues] = eig(v);
[mineigv,minindex] = min(diag(eigenValues));
x = eigenVectors(:,minindex);