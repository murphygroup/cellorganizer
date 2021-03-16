% 2011-12-08 tebuck: added "if isempty(maximum_area)" in case
% nuclear image is (intentionally) blank.

if isempty(nucstruct.channel_thr)

  if optimize && nuclear_channel_blank
  else

    nucstruct.channel_thr = graythresh( ...
        nucstruct.channel(nucstruct.channel>min(nucstruct.channel(:))))*127;
        
    nucstruct.channel_fg = nucstruct.channel;
    nucstruct.channel_fg(nucstruct.channel<=10) = 0;
%     nucstruct.channel_fg(nucstruct.channel<=nucstruct.channel_thr) = 0;
    
    tmp_bwl = bwlabel(nucstruct.channel_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>max(tmp_stats)/10);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    nucstruct.channel_fg(~tmp_bw) = 0;
    [maximum_area, maximum_index] = max(tmp_stats);
% $$$     maximum_area, maximum_index
% $$$     prctile(tmp_bwl(:), [0, 100])
    if ~isempty(maximum_area)
      nucstruct.channel_largest_object = tmp_bwl == maximum_index;
      nucstruct.channel_objectsizes = tmp_stats(tmp_idx);
    else
      nucstruct.channel_largest_object = tmp_bwl == 0;
      nucstruct.channel_objectsizes = 0;
    end
  
  end
  
end

