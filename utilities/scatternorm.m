function [ fh ] = scatternorm( data, class, plegend, ptitle, pxlabel, pylabel, pzlabel )
%Plots a scatter of 2 variables with elipses for gaussians fit to the
%emperical distribution
%
%var1 and var2 are cell arrays of equal length

if ~exist('class', 'var')
    class = ones(size(data,1),1);
end

fontsize = 16;

fh = figure;
% fh = figure('Units', 'Normalized', 'Position', [0.05, 0.05, 0.8, 0.8]);
set(gcf, 'Color', 'w');

hold on;



utype = unique(class);
colors = hsv(length(utype))* 0.8;

for i = 1:length(utype)
    tinds = class == utype(i);
    dat = data(tinds,:);
    
    color = repmat(colors(i,:), [length(dat), 1]);
    
    m = mean([dat]);
    c = cov([dat]);
    
    
    h(i) = plot_gaussian_ellipsoid(m, c);
    set(h(i), 'linewidth', 2)
    
    
    if size(dat, 2) == 2
        scatter(dat(:,1), dat(:,2), 300, color, '.')
        set(h(i), 'color', colors(i,:))
    else
        b = scatter3(dat(:,1), dat(:,2), dat(:,3), 300, color, '.');
        set(b, 'MarkerFaceColor', colors(i,:));
        
        set(h(i),'facecolor',colors(i,:))
        set(h(i),'facealpha', 0.25);
        set(h(i),'edgecolor','none');      
    end
    
    
    
end


if exist('plegend', 'var') & ~isempty(plegend)
    legend(h, plegend, 'location', 'NorthWest', 'FontSize', fontsize);
end


if exist('ptitle', 'var') & ~isempty(ptitle)
    title(ptitle, 'FontSize', fontsize);
end

if exist('pxlabel', 'var') & ~isempty(pxlabel)
    xlabel(pxlabel, 'FontSize', fontsize);
end

if exist('pylabel', 'var') & ~isempty(pylabel)
    ylabel(pylabel, 'FontSize', fontsize);
end

if exist('pzlabel', 'var') & ~isempty(pzlabel)
    zlabel(pzlabel, 'FontSize', fontsize);
end
set(gca,'FontSize',16)

end

