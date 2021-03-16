
feats = zeros(1, 1);
feats(1) = sum(protstruct.channel(:));
if optimize && protein_channel_blank
  feats(:) = nan;
end
names = {...
  'protTotalIntensity:total intensity'...
        };

slfnames = repmat({''}, 1, length(feats));
