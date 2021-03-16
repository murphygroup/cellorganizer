% Jieyue's previous code in commonScriptCalculateSet checked for 14
% features, i.e., inclusion of nuclear features, so nan if either
% channel is blank:
[default_names, default_feats, default_slfnames] = ml_imgfeatures(1, 1);
default_feats(:) = nan;
if optimize && (protein_channel_blank || nuclear_channel_blank)
% $$$   [names, feats, slfnames] = ml_imgfeatures(0, 0);
% $$$   [names, feats, slfnames] = ml_imgfeatures(1, 1);
% $$$   feats(:) = nan;
  names = default_names;
  feats = default_feats;
  slfnames = default_slfnames;
else
  [names, feats, slfnames] = ml_imgfeatures(protein_image, nuclear_image);
% $$$   if length(feats) ~= length(default_feats)
% $$$     feats = default_feats;
% $$$   end
end
