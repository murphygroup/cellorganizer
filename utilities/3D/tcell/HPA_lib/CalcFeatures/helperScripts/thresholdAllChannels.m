
thresholdProteinChannel

thresholdNucleusChannel

if isempty(tubstruct.channel_thr)

  if optimize && tubulin_channel_blank
  else

    tubstruct.channel_thr = graythresh( ...
        tubstruct.channel(tubstruct.channel>min(tubstruct.channel(:))))*255;
        
    tubstruct.channel_fg = tubstruct.channel;
    tubstruct.channel_fg(tubstruct.channel<=tubstruct.channel_thr) = 0;
    
    tmp_bwl = bwlabel(tubstruct.channel_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>5);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    tubstruct.channel_fg(~tmp_bw) = 0;
    tubstruct.channel_objectsizes = tmp_stats(tmp_idx);

  end
  
end

if isempty(erstruct.channel_thr)

  if optimize && er_channel_blank
  else

    erstruct.channel_thr = graythresh( ...
        erstruct.channel(erstruct.channel>min(erstruct.channel(:))))*255;
        
    erstruct.channel_fg = erstruct.channel;
    erstruct.channel_fg(erstruct.channel<=erstruct.channel_thr) = 0;
    
    tmp_bwl = bwlabel(erstruct.channel_fg>0, 4);
    
    tmp_props = regionprops(tmp_bwl,'Area');
    tmp_stats = zeros(length(tmp_props),1);
    for j = 1:length(tmp_props)
        tmp_stats(j,:) = tmp_props(j).Area;
    end
    tmp_idx = find(tmp_stats>5);
    
    tmp_bw = ismember(tmp_bwl, tmp_idx);
    erstruct.channel_fg(~tmp_bw) = 0;
    erstruct.channel_objectsizes = tmp_stats(tmp_idx);
    
  end
  
end

