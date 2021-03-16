
if optimize && er_channel_blank
else

    tmp_size = size(erstruct.channel);

    if isempty(erstruct.downsampled2x)
        erstruct.downsampled2x = imresize(erstruct.channel,[tmp_size(1) tmp_size(2)]/DSF);
    end

end
