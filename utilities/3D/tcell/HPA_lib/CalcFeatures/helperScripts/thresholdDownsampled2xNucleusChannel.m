
if isempty(nucstruct.downsampled2x_thr)

  if optimize && nuclear_channel_blank
  else

    nucstruct.downsampled2x_thr = graythresh( ...
        nucstruct.downsampled2x(nucstruct.downsampled2x>min(nucstruct.downsampled2x(:))))*127;
        
    nucstruct.downsampled2x_fg = nucstruct.downsampled2x;
    nucstruct.downsampled2x_fg(nucstruct.downsampled2x<=10) = 0;
%     nucstruct.downsampled2x_fg(nucstruct.downsampled2x<=nucstruct.downsampled2x_thr) = 0;
    
    tmp_bwl = bwlabel(nucstruct.downsampled2x_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>max(tmp_stats)/10);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    nucstruct.downsampled2x_fg(~tmp_bw) = 0;
    nucstruct.downsampled2x_objectsizes = tmp_stats(tmp_idx);

  end
  
end
