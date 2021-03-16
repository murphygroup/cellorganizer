function y = tz_evalfun(x,f,param)
%TZ_EVALFUN Evaluate a [parametric function].
%   Y = TZ_EVALFUN(X,F) returns the mapped values of X through the
%   [general function] F. The type of X depends on the function. If it
%   takes a vector or scalar value as a variable, then each row of X is a
%   sample of the variable. 
%   
%   Y = TZ_EVALFUN(X,F,PARAM) provides a way for function operation.
%   Currently scaling and translation are supported and they are only for
%   numeric fucntions. PARAM is a structure with the following two fields:
%       'xscale' - scale of the variables
%       'yscale' - scale of the returned values
%       'xtransl' - translation of the variables
%       'ytransl' - translation of the returned values
%
%       Notice: the scaling is done before translation
%
%   See also

%   10-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

% succ = 1;

warning(tz_genmsg('of','tz_evalfun','ml_evalfun'));

if ~exist('param','var')
    param = struct('init','on'); %set up an arbitrary field
end

if isfield(f,'transform')
    if isfield(f.transform,'xscale')
        if isfield(param,'xscale')
            param.xscale = param.xscale.*f.transform.xscale;
        else
            param.xscale = f.transform.xscale;
        end 
    end
    
    if isfield(f.transform,'yscale')
        if isfield(param,'yscale')
            param.yscale = param.yscale.*f.transform.yscale;
        else
            param.yscale = f.transform.yscale;
        end 
    end
    
    if isfield(f.transform,'xtransl')
        if isfield(param,'xtransl')
            param.xtransl = param.xtransl+f.transform.xtransl;
        else
            param.xtransl = f.transform.xtransl;
        end 
    end
    
    if isfield(f.transform,'ytransl')
        if isfield(param,'ytransl')
            param.ytransl = param.ytransl+f.transform.ytransl;
        else
            param.ytransl = f.transform.ytransl;
        end 
    end 
    
end
    
if isfield(param,'xscale')
    x = ml_multrow(x,param.xscale);
end

if isfield(param,'xtransl')
    x = ml_addrow(x,param.xtransl);
end

if isfield(f,'param')
    if iscell(f.param) & strcmp(f.funname,'ml_comfun')==0
        cmd = [f.funname '(x,f.param{:})'];
    else
        cmd = [f.funname '(x,f.param)'];
    end
else
    cmd = [f.funname '(x)'];
end

y = eval(cmd);

if isfield(param,'yscale')
    y = ml_multrow(y,param.yscale);
end

if isfield(param,'ytransl')
    y = ml_addrow(y,param.ytransl);
end

% 
% if succ==0
%     error('Function evaluation failed. Please check the input');
% end

