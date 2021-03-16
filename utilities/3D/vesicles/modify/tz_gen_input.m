function data =...
    tz_gen_input(T,rows,cols,errtol,feather,mode,pixel_delta_max, ...
    patch_delta_max,psr,psc,metric,mask)
%TZ_GEN_INPUT Modified from GEN_INPUT
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
% GEN_INPUT generates input data to hybridsynthesize.m
%
%   data = gen_input(T,rows,cols,errtol,feather,mode,
%                    pixel_delta_max,patch_delta_max,psr,psc,metric)
%
%   INPUT:
%   T - the input texture (upper left is at coord (1,1) in matlab)
%   rows/cols - the number of pixel rows and columns in the final image
%   errtol - added to miniumum error. assures diversity at the cost of 
%            possible bad matches. set to 0 to choose only best match, 1 to pick random
%   feather - set to 1 to smoothly blend the valid overlap pixels.
%             setting to 0 results in no feathering
%   mode - possible values are 'wrap' and 'nowrap' (defines the wrapping properties
%          of the input texture T). can cause excessive verbatim copying.
%   pixel_delta_max - error for which we validate pixels in the overlap region,
%                     setting to 0 resynthesizes the entire overlap region on a per 
%                     fragment basis. increasing pixel_delta_max, in the interval 
%                     of [0,1], trades off overlap feathering (see 'feather') for 
%                     per-pixel resynthesis. setting to 1 results in pure overlap 
%                     feathering.
%   patch_delta_max - maximal allowed error for patch lapping in [0,1]. 
%                     set to 1 for no adaptive splitting.
%   psr/psc - initial patch size (ps) for rows (r) and cols (c).
%   metric - error metric to be used for patch placement. possible values are
%            'simple', 'src', 'dst', 'sum', 'sqdiff'. if the string is unknown,
%            hybridsynthesize defaults to 'simple'
%
%   OUTPUT:
%   data - a structure holding all information to be processed by hybridsynthesize
%          this holds the above information and additionally:
%
%     1. all required fourier transforms (for the error metrics). combinations of
%        T, T^2, S and S^2
%     2. size of T (rows_t, cols_t)
%     3. an index_map for easy access of flattened (row major) indices.
%     4. the minimal patch size for adaptive splitting (min_patchsize)
%     5. the pixel overlap for per-pixel overlap resynthesis, both for k-coherence
%        (knn_overlap) and exhaustive (exh_overlap) search. note that the box-shaped
%        neighborhood around the target pixel has size [2*overlap+1]^2
%     6. the value for k, used in bulding the candidate set for k-coherence search.
%        smaller = faster, but possibly less accurate. set to k=1 for natural
%        textures in ashikhmin's sense
%     7. activation of overlap-resynthesis visualisation. set display_resynthesis
%        to 0 to disable
%     8. set how verbatim the algorithm should be. set to 0 to suppress all output
%
%   the returned datastructure is passed into hybridsynthesize.m
%   all values can (and should) be modified to allow the intended trade-off's.
%   see paper/thesis:
%
%   'Hybrid Texture Synthesis', Nealen A., Alexa M.
%   eurographics symposium on rendering 2003
%
%   'Hybrid Texture Synthesis', Nealen A.
%   Diplomarbeit (MSc thesis), Technische Universitaet Darmstadt
%

% source importance map (S)
S = frequencymap(T);
% alternatively uniform weighting
% S = ones(data.rows_t,data.cols_t);

T = im2double(T);    % convert input texture to double for arithmetic
S = im2double(S);    % convert importance map to double for arithmetic

% we generate an rgb output image, even if the input is grayscale (this makes special cases
% unecessary at a small cost). convert possible grayscale input texture to RGB color mode
T = gray2rgb(T);
S = gray2rgb(S);
data.T = T;

% most fourier transforms must only be computed at init time and are therefore stored
% in the input datastructure
ST     = S.*T;
Ts     = T.^2;
STs    = S.*Ts;
Ss     = S.^2;
TSs    = T.*Ss;
TsSs   = Ts.*Ss;

data.fftT    = fft2(T);         % fft of T
data.fftTs   = fft2(Ts);        % fft of T^2
data.fftS    = fft2(S);         % fft of S
data.fftST   = fft2(ST);        % fft of ST
data.fftSTs  = fft2(STs);       % fft of ST^2
data.fftSs   = fft2(Ss);        % fft of S^2
data.fftTSs  = fft2(TSs);       % fft of TS^2
data.fftTsSs = fft2(TsSs);      % fft of T^2*S^2

% build the input datastructure. note that none of this data is ever written to in
% the algorithm
data.rows_t          = size(T,1); 
data.cols_t          = size(T,2);
data.rows            = rows;
data.cols            = cols;
data.errtol          = errtol;
data.feather         = feather;
data.mode            = mode;
data.pixel_delta_max = pixel_delta_max;
data.patch_delta_max = patch_delta_max;
data.psr             = psr;
data.psc             = psc;
data.metric          = metric;

% index map. holds flattened indices for each position (r,c) in T.
data.index_map = zeros(data.rows_t,data.cols_t);
for r=1:data.rows_t,
    for c=1:data.cols_t,
        data.index_map(r,c) = (r-1)*data.cols_t + c;
    end
end

% some default values, which can be adjusted

% default box-shaped neighborhoods for exhaustive and k-coherence search in the 
% per-pixel resynthesis stage. the used neighborhood is (2*overlap + 1) wide
data.exh_overlap = 3;
data.knn_overlap = 2;

% the 'k' for k-coherence search. if this value is larger than the number of nearest 
% neighbors in the used knn structure (see build_knn.m) the algo defaults to the 
% available number of neighbors
data.k = 5;

% the smallest patch size allowed after splitting
data.min_patchsize = 8;

% display intermediate overlap-resynthesis results. setting to 0 disables this
data.display_resynthesis = 1;
data.verbatim = 0;

if exist('mask','var')
    data.mask = mask;
end

% perform some logic tests and give out warning(s)
num_row_patches = data.rows/data.psr;
num_col_patches = data.cols/data.psc;
if((num_row_patches-fix(num_row_patches)) ~= 0 |...
        (num_col_patches-fix(num_col_patches)) ~= 0),
    disp('WARNING: check that psr/psc is integer fraction of rows/cols.');
end

ps_row = psr;
while (ps_row > 2),
    ps_row = ps_row/2;
    if ((ps_row-fix(ps_row)) ~= 0),
        disp('WARNING: check that psr is a power of 2.');
        break;
    end
end
if (ps_row <= 2 & ps_row ~= 2 & ps_row ~= 1), 
    disp('WARNING: check that psr is a power of 2.'); 
end

ps_col = psc;
while (ps_col > 2),
    ps_col = ps_col/2;
    if ((ps_col-fix(ps_col)) ~= 0),
        disp('WARNING: check that psc is a power of 2.');
        break;
    end
end
if (ps_col <= 2 & ps_col ~= 2 & ps_col ~= 1), 
    disp('WARNING: check that psc is a power of 2.'); 
end
