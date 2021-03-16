function [ model ] = param2microtubule_model( cellparam, options )
%PARAM2MICROTUBULE_MODEL Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(cellparam)
    if ischar(cellparam{i})
        cellparam{i} = load(cellparam{i}, 'prot');
    end
end

if iscell(cellparam(1))
    cellparam = [cellparam{:}];
end
    
prot = [cellparam.prot];

mt_params = [vertcat(prot.n), vertcat(prot.mu_len), vertcat(prot.colli_min_number)];


model.type = 'network';
model.class = 'microtubule';
model.resolution = options.model.resolution;
model.version = '1.2';
model.parameters.mean = mean(mt_params);
model.parameters.cov = cov(mt_params);


end

