function [ nuclearparams, cellparams, proteinparams ] = getModelReportParams( model, param )
%GETMODELREPORTPARAMS 
%   Returns the names of the appropriate fields for each type of model

if ~exist('param', 'var') | isempty(param)
    param.includenuclear = 1;
    param.includecell = 1;
    param.includeprot = 1;
end

is3D = strcmpi(model.dimensionality, '3D');

proteinparams = {};
nuclearparams = {};
cellparams = {};

if is3D
    if param.includeprot
        if strcmpi(model.proteinModel.type,'gmm')
            proteinparams = {'size.x.mu', 'size.x.sigma', ...
                'size.y_x.a1', 'size.y_x.b1', 'size.y_x.a2', 'size.y_x.b2', ...
                'size.z_x.a1', 'size.z_x.b1', 'size.z_x.a2', 'size.z_x.b2', ...
                'frequency.mu', 'frequency.sigma' ...
            %    , 'position.beta(1)', 'position.beta(2)', 'position.beta(3)', ...
            %    'position.beta(4)', 'position.beta(5)', 'position.beta(6)'
                };
        end
        if strcmpi(model.proteinModel.type,'spharm_obj')
            %proteinparams={'X','coeff','train_score','all_spharm_descriptors',...
            %    'ppm.full.covars_max','ppm.full.covars_min','ppm.full.param'};
        end
    end
    
    if param.includenuclear
        nuclearparams = {'height.stat.mu', 'height.stat.sigma'}; 
    end
    
else
    
    if param.includeprot
        proteinparams = {'objectModel.stat.beta', ...
            'objectModel.intensStatModel.mu', ...
            'objectModel.intensStatModel.sigma', ...
            'objectModel.numStatModel.alpha', ...
            'objectModel.numStatModel.beta'};
    end
    
    if param.includenuclear
        nuclearparams = {}; 
    end
end

%this is not set
cellparams = {};

end

