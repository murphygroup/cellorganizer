function tz_plotfun(f,param)
%TZ_PLOTFUN Plot a function.
%   TZ_PLOTFUN(F,PARAM) plots the [parametric function] F according to the
%   structured parameter PARAM, which has the following fields (* means the 
%   filed is necessary to be initialized):
%       * range: a kx2 vector, where k is the number of variables.
%       npoint: a kx1 vector for number of points for each variable
%       evalfun: the third parameter for tz_evalfun
%       pol2cart: transform polar to Cartesian coordinates
%       title: title of the plot
%       'plot' - a cell array secifying plot options (default
%           {})
%   Notice: this function only supports 1D or 2D function.
%   
%   See also

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ~isfield(param,'range')
    error('You must specify the range of the function');
end

k = size(param.range,1);

if k>2
    error('This function only supports 1D or 2D function.');
end

funname = strrep(f.funname,'_','\_');
param = ml_initparam(param,struct('npoint',repmat(100,k,1), ...
    'evalfun',[],'pol2cart',0,'title',funname,'plot',{{}}));

if k==1
    xs = tz_ppoints(param.range(1),param.range(2),param.npoint);
    plot(xs,ml_evalfun(xs',f,param.evalfun),param.plot{:});

else
    xs = tz_ppoints(param.range(1,1),param.range(1,2),param.npoint);
    ys = tz_ppoints(param.range(2,1),param.range(2,2),param.npoint);
    [xs,ys]=meshgrid(xs,ys);
    gridSize = size(xs);
    
    z = ml_evalfun([xs(:),ys(:)],f,param.evalfun);
    
    
    if param.pol2cart==1
        [xs,ys,z] = pol2cart(ys(:),xs(:),z(:));
        xs = reshape(xs,gridSize);
        ys = reshape(ys,gridSize);
    end
    z = reshape(z,gridSize);
    surf(xs,ys,z);
    zlabel('z');
end

xlabel('x');
ylabel('y');
if ~isempty(param.title)
    title(param.title);
end

