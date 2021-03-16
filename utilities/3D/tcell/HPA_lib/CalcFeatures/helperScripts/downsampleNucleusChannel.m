
if optimize && nuclear_channel_blank
else

    tmp_size = size(nucstruct.channel);

    if isempty(nucstruct.downsampled2x)
        nucstruct.downsampled2x = imresize(nucstruct.channel,[tmp_size(1) tmp_size(2)]/DSF);
    end

end
