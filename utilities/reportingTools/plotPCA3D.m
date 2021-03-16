function [ cumLatent, score ] = plotPCA3D(data, labels, project, cm_handle)
%Plots the first 3 prinicple componenets

    
    if ~exist('project', 'var') | isempty(project)
        project = false;
    end
    
    if ~exist('labels', 'var')
        labels = 1:size(data,1);
    end
    
    rminds = any(isnan(data),2);
    
    data(rminds,:) = [];
    labels(rminds) = [];
    
    [ulabels,~,labelinds] = unique(labels);
    
    if ~exist('cm_handle', 'var')
        colors = jet(length(ulabels))* 0.8;
    else
        colors = cm_handle(length(ulabels)).*0.8;
    end
    
    data = zscore(data);
    
    [coeff, score, latent] = pca(data, 'Centered', false);
    cumLatent = cumsum(latent./sum(latent));

    %find the PCA'd dimensions that account for 0.95% of the variance
    inds = find(cumLatent >= 0.95);
    
    if size(score,2) < 3;
        score = padarray(score, [0,3-size(score,2)],'post');
    end
    
    hold on
    scatter3(score(:,1),score(:,2), score(:,3),25, colors(labelinds,:), 'filled')
    view(3)
    axis tight
    a = axis;
    
    if exist('project', 'var') && project
        for i = 1:size(score,1)
           line([a(2), score(i,1)]',[score(i,2) score(i,2)]', [score(i,3) score(i,3)]', 'Color', [0.75, 0.75, 0.75]);
           line([score(i,1), score(i,1)]',[a(4) score(i,2)]', [score(i,3) score(i,3)]', 'Color', [0.75, 0.75, 0.75]);
           line([score(i,1), score(i,1)]',[score(i,2) score(i,2)]', [a(5) score(i,3)]', 'Color', [0.75, 0.75, 0.75]);
        end
        scatter3(score(:,1),score(:,2), score(:,3), [], colors(labelinds,:),'.')
    else
        grid on
    end
    hold off
    
    
    xlabel('pc1', 'FontSize', 16);
    ylabel('pc2', 'FontSize', 16);
    zlabel('pc3', 'FontSize', 16);
    
%     grid on

%     plegend = cellstr(num2str(ulabels'));
%     legend(a, plegend, 'location', 'NorthWest', 'FontSize', 16);
%     legend off

    set(gcf, 'Color', 'w')
    
end