function plotPCAGauss(feats, times, title)

if ~exist('title', 'var') | isempty(title)
    title = '';
end

    [~, score] = princomp(feats);

    utimes = unique(times);

    scatternorm(score(:,1:3), times, cellstr(num2str(utimes')), title, 'pc1', 'pc2', 'pc3')
end