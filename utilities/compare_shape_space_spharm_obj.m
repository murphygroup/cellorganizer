function [ consolidated_model ] = compare_shape_space_spharm_obj( models,fileID,options )

% used for comparisons of both cell/nuclear shape models and for spharm_obj models

% 8/13/2022 R.F.Murphy add comparison of paired models
% 3/21/2023 R.F.Murphy fix thumbnail saving

    consolidated_model=consolidateshapespacemodels(models{1},models{2});

% if the samples in the two models are paired (i.e., different models run on the
% same set of images), analyze their pairwise distances in the reduced space
    if options.paired
        size1=size(models{1}.X,1);
        size2=size(models{2}.X,1);
        if size1~=size2
            warning('Unable to compare samples pairwise: size mismatch');
        else
            for isamp = 1:size1
                scorediff(isamp) = pdist2(consolidated_model.train_score(isamp,:),consolidated_model.train_score(size1+isamp,:));
            end
            header2html(fileID,'Distances between paired shapes');
            text2html(fileID,['Mean=',num2str(mean(scorediff))]);
            text2html(fileID,['Min=',num2str(min(scorediff))]);
            text2html(fileID,['Max=',num2str(max(scorediff))]);
            consolidated_model.paired_score_diff = scorediff;
        end
    end

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
    fP = f.Position;
    f.Position = [fP(1) fP(2) round(fP(3)/2) round(fP(4)/2)];
    fname2 = sprintf('show_shape_space_%d_thumbnail.png', i);
    saveas( f, fname2, 'png' ); 3/21/2023
    img2html(fileID,fname,fname2,'Joint shape space from SPHARM-RPDM models.');
end
end