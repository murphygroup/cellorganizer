
function [ consolidated_model ] = compare_shape_space_spharm_obj( models,fileID,options )
consolidateX = [];
idx = [0];
for i=1:length(models)
    X1=models{i}.X;
    consolidateX = vertcat(consolidateX, X1);
    idx = [idx idx(end)+size(models{i}.X, 1)];
end

%default value in cellorganizer
latent_dim=100;

[scales,coeff,score,latent,tsquared,explained,mu,train_score_consolidate,train_explained,train_coeff_consolidate] = CalculatePCA(consolidateX, latent_dim);

for i=1:length(models)-1
    model1 = models{i};
    model2 = models{end};
    spharm1=model1.all_spharm_descriptors;
    spharm2=model2.all_spharm_descriptors;

    %Create new model with required features 
%     consolidated_model.coeff=vertcat(train_coeff_consolidate(idx(i)+1:idx(i+1),:), train_coeff_consolidate(idx(end-1)+1:idx(end),:));
    consolidated_model.X=consolidateX;
    consolidated_model.train_score=vertcat(train_score_consolidate(idx(i)+1:idx(i+1),:), train_score_consolidate(idx(end-1)+1:idx(end),:));
    consolidated_model.numimgs1=model1.numimgs;
    consolidated_model.numimgs2=model2.numimgs;
    consolidated_model.numimgs=model1.numimgs+model2.numimgs;
    consolidated_model.all_spharm_descriptors=cat(3,spharm1,spharm2);
    if isfield(model1, 'hausdorff_distances') && isfield(model2, 'hausdorff_distances')
        consolidated_model.hausdorff_distances=horzcat(model1.hausdorff_distances,model2.hausdorff_distances);
    end
    if ismember(model1.components,model2.components)
        consolidated_model.components=model2.components;
    else
        consolidated_model.components={};
        consolidated_model.components={model1.components,model2.components};
    end
    f = figure('visible','off');
    updateLabels(consolidated_model,options);
    fname = sprintf('show_shape_space_%d.png', i);
    saveas( f, fname, 'png');
    I = imread(fname);
    I = imresize( I, 0.50 );
    fname2 = sprintf('show_shape_space_%d_thumbnail.png', i);
    imwrite( I, fname2 );
    img2html(fileID,fname,fname2,'Joint shape space from SPRM-RPDM models.');
end
end

