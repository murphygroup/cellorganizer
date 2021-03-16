
if isempty(protstruct.downsampled2x_mthr)
    protstruct.downsampled2x_thr = graythresh( ...
        protstruct.downsampled2x(protstruct.downsampled2x>min(protstruct.downsampled2x(:))))*255;

    protstruct.downsampled2x_mthr = graythresh( ...
        protstruct.downsampled2x(protstruct.downsampled2x>protstruct.downsampled2x_thr))*255;
        
    protstruct.downsampled2x_fg = protstruct.downsampled2x;
    protstruct.downsampled2x_fg(protstruct.downsampled2x<=protstruct.downsampled2x_mthr) = 0;
    
    tmp_bwl = bwlabel(protstruct.downsampled2x_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>max(tmp_stats)/10);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    protstruct.downsampled2x_fg(~tmp_bw) = 0;
    protstruct.downsampled2x_objectsizes = tmp_stats(tmp_idx);
    
    tmp_idx = find(tmp_stats>max(tmp_stats(:))/5);
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    tmp_bw = imdilate( tmp_bw, strel('disk',4));
    protstruct.downsampled2x_large_fg = tmp_bw;
    
    tmp_bwl = bwlabel(tmp_bw, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end

    protstruct.downsampled2x_large_objectsizes = tmp_stats;
end
