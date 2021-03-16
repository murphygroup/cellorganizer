function [lambda,ordering] = shape_interp_track(x,c)
%SHAPE_INTERP_TRACK calculates the interpolation order to generate a shape
%inside the high-dimensional triangle using the vertices
%
%Adapted from Peng et al 2009.
%http://murphylab.web.cmu.edu/publications/158-peng2009.pdf, see Methods 2.4

% 2013-02-11 tebuck: added comments to clarify.



[m,n] = size(x);
if m ~= n-1 || n < 3
    % error('Not a triangle!')
    error('x is %d x %d but should be (n-1) x n to be a simplex!', m, n)
end

lambda = zeros(1,n-1);
ordering = zeros(1,n);
norder = 1:n;
for i = 1:n-2
    % Order the vertices by distance from the query point c:
    distseq = sum((x-repmat(c,1,size(x,2))).^2,1);
    % Pick the vertex at minimum distance for this interpolation:
    [mindist,minidx] = min(distseq);
    % x_(k+1) from the paper:
    x1 = x(:,minidx);
    x(:,minidx) = [];
    ordering(i) = norder(minidx);
    norder(minidx) = [];
    % Calculate segmentation ratio (lambda_k apparently):
    Num = x - repmat(x(:,1),1,size(x,1));
    Num(:,1) = c - x(:,1);
    Den = Num;
    Den(:,1) = c - x1;
    lambda(i) = - det(Num) / det(Den);
    % Dimension reduction (by projection onto a face of the current simplex):
    w = [c+lambda(i)*(c-x1),x]';
    [coef,w] = pca(w);
    w(:,end) = [];
    % c_(k-1) from the paper:
    c = w(1,:)';
    x = w(2:end,:)';
end

dist = abs(x-c);

[mindist,minidx] = min(dist);
dist(minidx) = [];
lambda(n-1) = dist / mindist;
ordering(n-1) = norder(minidx);
norder(minidx) = [];
ordering(n) = norder;
