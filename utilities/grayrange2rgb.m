function [ colors_out ] = grayrange2rgb( grayrange_all, gmin, gmax, cm )
%takes a range of gray values and returns those positions on a colormap

npts = length(grayrange_all);

validinds = ~isnan(grayrange_all) & ~isinf(grayrange_all);

grayrange = grayrange_all(validinds);


if ~exist('gmin', 'var') | isempty(gmin)
    gmin = min(grayrange);
end

if ~exist('gmax', 'var') | isempty(gmax)
    gmax = max(grayrange);
end

if ~exist('cm', 'var') | isempty(cm)
    cm = @jet;
end



bit_depth = 2^16;

colorinds = ((grayrange - gmin)/(gmax-gmin)) * (bit_depth-1) + 1;

%interpolate between 16-bit colors
colors = cm(bit_depth);


colorinds(colorinds < 1) = 1;
colorinds(colorinds > bit_depth) = bit_depth;

colors_out = zeros(npts, 3);

colors_out(validinds,:) = colors(round(colorinds),:);



end

