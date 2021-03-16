
thresholdDownsampled2xProteinChannel

thresholdDownsampled2xNucleusChannel

if isempty(tubstruct.downsampled2x_thr)

  if optimize && tubulin_channel_blank
  else

    tubstruct.downsampled2x_thr = graythresh( ...
        tubstruct.downsampled2x(tubstruct.downsampled2x>min(tubstruct.downsampled2x(:))))*255;
        
    tubstruct.downsampled2x_fg = tubstruct.downsampled2x;
    tubstruct.downsampled2x_fg(tubstruct.downsampled2x<=tubstruct.downsampled2x_thr) = 0;
    
    tmp_bwl = bwlabel(tubstruct.downsampled2x_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>max(tmp_stats)/10);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    tubstruct.downsampled2x_fg(~tmp_bw) = 0;
    tubstruct.downsampled2x_objectsizes = tmp_stats(tmp_idx);

  end
  
end

if isempty(erstruct.downsampled2x_thr)

  if optimize && er_channel_blank
  else

    erstruct.downsampled2x_thr = graythresh( ...
        erstruct.downsampled2x(erstruct.downsampled2x>min(erstruct.downsampled2x(:))))*255;
        
    erstruct.downsampled2x_fg = erstruct.downsampled2x;
    erstruct.downsampled2x_fg(erstruct.downsampled2x<=erstruct.downsampled2x_thr) = 0;
    
    tmp_bwl = bwlabel(erstruct.downsampled2x_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>max(tmp_stats)/10);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    erstruct.downsampled2x_fg(~tmp_bw) = 0;
    erstruct.downsampled2x_objectsizes = tmp_stats(tmp_idx);

  end
  
end
