function [index,distance] = build_knn(T,ov,k,retain,mode)
%
% |----------------------------------------------------------|
% | Hybrid Texture Synthesis MATLAB package                  |
% |                                                          |
% | Author: Andrew Nealen                                    |
% |         Discrete Geometric Modeling Group                |
% |         Technische Universitaet Darmstadt, Germany       |
% |                                                          |
% | Note:   This is part of the prototype implementation     |
% |         accompanying our paper/my thesis                 |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen and M. Alexa |
% |         Eurographics Symposium on Rendering 2003         |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen              |
% |         Diplomarbeit (MSc), TU Darmstadt, 2003           |
% |                                                          |
% |         See the paper/thesis for further details.        |
% |----------------------------------------------------------|
%
% BUILD_KNN build k-nearest neighbors for all pixels of input texture T
%
%   [index, distance] = build_knn(T,ov,k,retain,mode)
%
%   NOTE: this function requires openTSTool's nn_search and nn_prepare, both
%         avilable under the GPL from http://www.physik3.gwdg.de/tstool/
%         where detailed installation instructions are available
%
%   input:
%   T      - the input texture of size RxC=N (boundaries wrap in 2D)
%   ov     - the pixel overlap (the nxn neighborhood for feature vector 
%            construction with n = 2*ov + 1)
%   k      - number of nearest neighbors to search for (and return)
%   retain - amount of variance to retain during PCA (in %) (optional)
%   mode   - set to 1 for grayscale reduction (optional)
%
%   output:
%   index    - an (Nxk) sized matrix with the k-nearest neighbor indices
%              for each sample in SCANLINE (row major) order
%   distance - an (Nxk) sized matrix with the corresponding euclidian
%              distances
%
%   NOTE: the result is a flattened version of the input 2D matrix with
%         SCANLINE ordering. to get the vector of (flat) nearest neighbor 
%         indices at pixel p = f(r,c) we access 'index(p, :)' with
%
%         p = f(r,c) = (r-1)*C+c
%
%         with C = 'num cols in T'.
%         these indices are flattened, so to get the actual (r,c) index
%         of the flattened index p we compute
%
%         r = floor((p-1)/C) + 1, and
%         c = p - (r-1)*C
%
%         with r = row index and c = col index
%         these transforms are commonly used in hybridsynthesize.m related code
%

T = im2double(T); % gives us values in [0,1]

% default to retain 97% of the variance during PCA
pca_retain = 97;
if (exist('retain')),
    pca_retain = retain;
end

% default to RGB mode, set to grayscale on user request
if (exist('mode') & mode == 1),
    if (size(T,3) > 1),
        T = rgb2gray(T);
    end
end

% color dimension of input image
col_dim = size(T,3);
col_dim_minus_one = col_dim-1;

% build the input to nn_prepare and nn_search
R = size(T,1); C = size(T,2);     % num rows (R) and cols (C) in T
N = R * C;                        % number of points in search space (N)
nbhd = 2*ov+1;                    % quadratic nbhd size
D = (nbhd^2) * col_dim;           % dimensionality (D) of each search vector
pointset = zeros(N, D);           % vector of length D for each pixel position 1...N

% ----------------------------------------------------------
% build pointset from T as input to nn_prepare and nn_search
% ----------------------------------------------------------

% shift it circulary, so entire nbhd is in upper left. also, the index of 
% the center pixel will always be held at (b_ind,b_ind), and we circularly
% shift T in each step (in scanline order)
T = circshift(T,[ov ov]);         
center_ind = 1+ov;

disp(strcat('building set of N=', num2str(N), ...
    ', feature vectors with dimension D=', num2str(D)));

tic

% iterate over all pixels in T in scanline order. 
% construct feature vectors ('pointset').
point_index = 0;
for r=1:R,
    for c=1:C,
        % current point index. row major (scanline order). runs from 1...N
        % point_index = (r-1)*C+c;
        point_index = point_index + 1;
        % for each vector V(point_index), flatten nbhd around current pixel
        pointset(point_index,1:D) = reshape(T(1:nbhd,1:nbhd,1:col_dim),1,D);
        % end of this column, circshift one column to the left, so current pixel evaluated
        % is one pixel to the right of previous evaluated pixel
        T = circshift(T,[0 -1]);
    end
    % end of this row, circshift one row up, so current row evaluated
    % is one pixel below the previous row evaluated
    T = circshift(T,[-1 0]);
end

toc

tic

% ----------------------------------------------------------
% perform PCA on pointset to reduce dimensionality
% comment this out if you 'really' want to work without PCA
% warning: no PCA = SLOW!
% ----------------------------------------------------------
disp('reducing feature vector dimensionality by PCA');
% retain 97% of the original variance as proposed by liang et al. in 
% 'real time texture synthesis by patch based sampling' paper.
[rlvm,frvals,frvecs,pointset] = tz_pca(pointset,'normalized',pca_retain,1);

% ----------------------------------------------------------
% prepare and search in pointset
% ----------------------------------------------------------
disp('preparing set of feature vectors for knn search');

% prepare pointset
atria = nn_prepare(pointset);

disp(strcat('computing knn with k=', num2str(k),', num feature vectors N=', num2str(N),...
    ', with dimension D=', num2str(D), ', reduced to D_pca=', num2str(size(pointset,2))));

% compute (exact) k-nearest neighbors
% the fifth paramtere to nn_search defines the accuracy of the search (0 = exact)
% see documentation accompanying nn_search in openTSTools to tweak this
[index,distance] = nn_search(pointset, atria, [1:N], k, 0);

disp('done');

toc