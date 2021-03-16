function spfeat = tp_spcylfeat(surfmap,param)
% TP_CYLFEAT extract b-spline surface coefficients as generative
% features of cylindrical form of 3D surface

if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,...
            struct('xorder',3,'yorder',4,...
            'xknots',2,'yknots',6,...
            'gapweight',1000) );
%         struct('xorder',3,'yorder',4,...
%         'xknots',[0 0 0 .5 1 1 1],...
%         'yknots',[0 0 0 0 1/4 3/8 1/2 5/8 3/4 1 1 1 1],...
%         'gapweight',1000));

[H,W] = size(surfmap);
x = 0:(1/(H-1)):1;
y = -pi:(2*pi/(W-1)):pi;
w = ones(size(y));
w([1,end]) = param.gapweight;


spfeat = spap2({param.xknots,param.yknots},...
    [param.xorder,param.yorder],{x,y},surfmap,{[],w});

spfeat.coefs = squeeze(spfeat.coefs);
spfeat.height = H-1;