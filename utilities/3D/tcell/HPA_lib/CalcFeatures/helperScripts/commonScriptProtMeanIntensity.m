
feats = zeros(1, 1);
if optimize && (protein_channel_blank || tubulin_channel_blank)
  feats(:) = nan;
else
  relevant_pixels = tubstruct.channel_fg > 0;
  feats(1) = sum(protstruct.channel(relevant_pixels)) ./ sum(relevant_pixels(:));
end
names = {...
  'protMeanIntensity:mean intensity'...
        };

slfnames = repmat({''}, 1, length(feats));
