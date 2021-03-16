function cellimg = tp_gencellshape(model,nuclei,param)
% TP_GENCELLSHAPE generates a cell shape from the eigen shape model of the
% radius ratio of the nuclei to the cell membrane.

if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

% param = ml_initparam(param,...
%     struct('xsize',1024,'ysize',1024,'samp_rate',360,'alpha',0.002));
param = ml_initparam(param,...
    struct('xsize',1024,'ysize',1024,'samp_rate',360,'alpha',0.002));

xcenter = param.xsize / 2;
ycenter = param.ysize / 2;

norm_std.name = 'mvn';
norm_std.mu = zeros(1,size(model.modeShape.const,2));
norm_std.sigma = eye(size(model.modeShape.const,2));
syn_score = ml_rnd(norm_std);
while tp_pvalue(syn_score,norm_std) < param.alpha
    syn_score = ml_rnd(norm_std);
end
syn_score = syn_score .* sqrt(model.eigens.stat');

shape_vec = model.meanShape.const + (model.modeShape.const * syn_score')';

height = length(shape_vec)/360;
nuc2cell_ratio = reshape(shape_vec,[height 360]);
%nuc2cell_ratio = reshape(shape_vec,[18 360]);

% Up-sampling
MEAN_HEIGHT = 85;
t = 0:(height-1)/MEAN_HEIGHT:height-1;
nuc2cell_ratio_interp = zeros(MEAN_HEIGHT+1,360);
for i = 1:360
    nuc2cell_ratio_interp(:,i) = interp1(0:height-1,nuc2cell_ratio(:,i),t);
end

nucimgsize = nuclei.nucimgsize;
% cellheight = nucimgsize(3);
nuc2cell_ratio_interp(:,end+1) = nuc2cell_ratio_interp(:,1);

nucsurf = nuclei.nucsurf;
delta = 2*pi/param.samp_rate;
Phi = -pi:delta:pi;

cellimg = zeros(nucimgsize(1),nucimgsize(2),MEAN_HEIGHT+1);

for i = 1:MEAN_HEIGHT+1
    [x,y] = pol2cart(Phi,nucsurf(i,:)./nuc2cell_ratio_interp(i,:));
    x = x + xcenter;
    y = y + ycenter;
    sliceimg = zeros(param.ysize,param.xsize);
    for t = 1:length(x)-1
        rpts = round(linspace(y(t),y(t+1),50));
        cpts = round(linspace(x(t),x(t+1),50));
        try
            index = sub2ind([param.ysize,param.xsize],rpts,cpts);
            sliceimg(index) = 255;
        catch
            warning(['CellOrganizer: Out of range subscript ' num2str(t) ' when adding pixel to image'])
        end
    end
    cellimg(:,:,i) = imfill(sliceimg,'holes');
end
