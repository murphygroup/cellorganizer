function [medfeat,heifeat] = tp_medplanefeat(medplane,height,mask,param)
% TP_MEDPLANEFEAT extract b-spline surface coefficients as generative
% features of medial plane and height

if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,...
    struct('xorder',4,'yorder',4,'xknots',1,'yknots',1));

[H,W] = size(medplane);
x = (0:H-1)/H;
y = (0:W-1)/W;

spmed = spap2({param.xknots,param.yknots},...
    [param.xorder,param.yorder],{x,y},medplane);
sphei = spap2({param.xknots,param.yknots},...
    [param.xorder,param.yorder],{x,y},height);

medfeat = struct('coefs',squeeze(spmed.coefs));
heifeat = struct('coefs',squeeze(sphei.coefs));