%if optimize && any(strcmpi(naming_convention.blank_channels, 'protein'))
if optimize && protein_channel_blank
  [names, feats, slfnames] = ml_tas(0, 0);
  %feats = nan(size(feats)); 
  feats(:) = nan;
else
  [names, feats, slfnames] = ml_tas(prot, 0);
end
