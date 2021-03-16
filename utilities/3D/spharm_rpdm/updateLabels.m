function [ labels ] = updatelabels( model, options)
%Addition : Hanxi Xiao and Serena Abraham
%Update Labels corresponding to each model
train_score = model.train_score;
nimgs = size(train_score, 1);


if isfield(model,'numimgs1') && isfield(model,'numimgs2')
       labels = ones(nimgs, 1);
       labels(model.numimgs1+1:end)=2*labels(model.numimgs1+1:end);
else
       labels = ones(nimgs, 1);
end
show_SPHARM_RPDM_Shape_Space_Figure(model,labels,options);

end
