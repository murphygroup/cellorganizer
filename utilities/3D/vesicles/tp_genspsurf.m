function instance = tp_genspsurf(model,param)
% Generate spline surface instance from a statistical spline surface model

% July 23, 2012 R.F. Murphy Add code for testing with mean shape
%                           (uncomment to use)
% June 29, 2013 R.F. Murphy Add param.generatemeanshape to get mean shape
% July 1, 2013 R.F. Murphy  Add default param.generatemeanshape=false

if nargin < 2
    param = [];
end    
param = ml_initparam(param,struct('alpha',0.002,'generatemeanshape',false));



% Side surface coefficients
f_coef = model.surface.stat;
if param.generatemeanshape
    height = round(model.height.stat.mu);
    
    coefs = f_coef.mu;
    tp_pvalue(coefs,f_coef);
else
    
    % Height of the cell
    f_height = model.height.stat;

    height = round(ml_rnd(f_height));
    while tp_pvalue(height,f_height) < param.alpha
        height = round(ml_rnd(f_height));
    end
    
    
    coefs = ml_rnd(f_coef);
    while tp_pvalue(coefs,f_coef) < param.alpha
        coefs = ml_rnd(f_coef);
    end
end

number = model.surface.number;
number(2) = number(2) - 1;
coefs = reshape(coefs,number);
coefs(:,end+1) = coefs(:,1);

% Formalize the generated shape parameters into spline structure
instance.coefs = coefs;
instance.form = model.surface.form;
instance.knots{1} = [zeros(1,model.surface.order(1)) ...
                model.surface.constknot_h ...
                ones(1,model.surface.order(1))];
instance.knots{2} = pi*[-1*ones(1,model.surface.order(2)) ...
                2*model.surface.constknot_phi-1 ...
                ones(1,model.surface.order(2))];
instance.number = model.surface.number;
instance.order = model.surface.order;
instance.dim = 1;
instance.height = height;
