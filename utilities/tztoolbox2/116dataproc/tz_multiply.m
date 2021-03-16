function x = tz_multiply(x1,x2)
%TZ_MULTIPLY N-D matrix multiplication
%   X = TZ_MULTIPLY(X1,X2) multiply two multiple dimentional matrix. The
%   returned value is a L*M*...*N*U*V*...*W matrix if X1 is a 
%   L*M*...*N*K matrix ans X2 is a K*U*V*...*W matrix.
%   
%   See also

%   03-Nov-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ndims(x1) ==1 | ndims(x2)==1 | (ndims(x1)<=2 & ndims(x2)<=2)
    x = x1*x2;
end

size1 = size(x1);
size2 = size(x2);

x1 = reshape(x1,[prod(size1(1:end-1)),size1(end)]);
x2 = reshape(x2,[size2(1),size2(2:end)]);

x = x1*x2;

x = reshape(x,[size1(1:end-1),size2(2:end)]);