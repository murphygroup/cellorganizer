
if optimize && tubulin_channel_blank
else

    tmp_size = size(tubstruct.channel);

    if isempty(tubstruct.downsampled2x)
        tubstruct.downsampled2x = imresize(tubstruct.channel,[tmp_size(1) tmp_size(2)]/DSF);
    end

end
