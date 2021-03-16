
if isempty(protstruct.channel_mthr)
    protstruct.channel_thr = graythresh( ...
        protstruct.channel(protstruct.channel>min(protstruct.channel(:))))*255;

    protstruct.channel_mthr = graythresh( ...
        protstruct.channel(protstruct.channel>protstruct.channel_thr))*255;

    protstruct.channel_fg = protstruct.channel;
    protstruct.channel_fg(protstruct.channel<=protstruct.channel_mthr) = 0;

    tmp_bwl = bwlabel(protstruct.channel_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end

    tmp_idx = find(tmp_stats>5);
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    protstruct.channel_fg(~tmp_bw) = 0;
    protstruct.channel_objectsizes = tmp_stats(tmp_idx);

    tmp_idx = find(tmp_stats>max(tmp_stats(:))/5);
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    tmp_bw = imdilate( tmp_bw, strel('disk',4));
    protstruct.channel_large_fg = tmp_bw;
    
    tmp_bwl = bwlabel(tmp_bw, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end

    protstruct.channel_large_objectsizes = tmp_stats;
end
