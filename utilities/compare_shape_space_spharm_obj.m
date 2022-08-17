function [ consolidated_model ] = compare_shape_space_spharm_obj( models,fileID,options )

% used for comparisons of both cell/nuclear shape models and for spharm_obj models

% 8/13/2022 R.F.Murphy add comparison of paired models

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

    f = figure('visible','off');
    updateLabels(consolidated_model,options);
    saveas( f, sprintf('show_shape_space.png', 'png'));
    I = imread( 'show_shape_space.png' );
    I = imresize( I, 0.50 );
    imwrite( I, 'show_shape_space_thumbnail.png' );
    img2html(fileID,'show_shape_space.png','show_shape_space_thumbnail.png','Joint shape space from SPRM-RPDM models.');
end