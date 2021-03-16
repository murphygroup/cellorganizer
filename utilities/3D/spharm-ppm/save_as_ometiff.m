function save_as_ometiff(cellimg,dnaimg,protimg,parameters,sfp) % remember to put mask path to parameters
    [rows, cols, ~] = size(cellimg);
    for z=1:1:size(cellimg,3)
        omeimg(:,:,z,1,1) = imresize(cellimg(:,:,z), [rows cols]);
    end
    for z=1:1:size(dnaimg,3)
        omeimg(:,:,z,2,1) = imresize(dnaimg(:,:,z), [rows cols]);
    end
    
    for z=1:1:size(protimg,3)
        omeimg(:,:,z,3,1) = imresize(protimg(:,:,z), [rows cols]);
    end
    metadata = createOMEXMLMetadata(omeimg,parameters,'XYZCT');
    %metadata.PhysicalSizeX
    if isfield( parameters, 'PhysicalSizeX' )
        pixelSize = ome.units.quantity.Length( ...
            java.lang.Double(parameters.PhysicalSizeX), ome.units.UNITS.MICROM);
        metadata.setPixelsPhysicalSizeX(pixelSize, 0);
    else
        warning('PhysicalSizeX not set. Exiting method');
        return
    end

    %metadata.PhysicalSizeY
    if isfield( parameters, 'PhysicalSizeY' )
        pixelSize = ome.units.quantity.Length( ...
            java.lang.Double(parameters.PhysicalSizeY), ome.units.UNITS.MICROM);
        metadata.setPixelsPhysicalSizeY(pixelSize, 0);
    else
        warning('PhysicalSizeX not set. Exiting method');
        return
    end

    %metadata.PhysicalSizeZ
    if isfield( parameters, 'PhysicalSizeZ' )
        pixelSize = ome.units.quantity.Length( ...
            java.lang.Double(parameters.PhysicalSizeZ), ome.units.UNITS.MICROM);
        metadata.setPixelsPhysicalSizeZ(pixelSize, 0);
    else
        disp('PhysicalSizeZ not set. Assuming you are attempting to create a 2D image');
    end

    %metadata.list_of_channel_labels
    if ~isfield( parameters, 'list_of_channel_labels')
        warning( 'Mandatory parameter list_of_channel_labels does not exist. Exiting method.' );
        return

%     elseif length(parameters.list_of_channel_labels) ~= length(list_of_input_images)
%         disp('auas: ');
%         disp(length(list_of_input_images));
%         disp(length(parameters.list_of_channel_labels));
%         warning('The number of input images must be equal to the number of labels. Exiting method.' );
%         return
%     else
    end

    %metadata.list_of_channel_labels
    for index=1:1:length(parameters.list_of_channel_labels)
        channel_name = parameters.list_of_channel_labels{index};
        channel_index = index-1;
        metadata.setChannelName( java.lang.String(channel_name), 0, channel_index )
    end
    disp(['Saving to ',sfp])
    bfsave( omeimg, sfp, 'metadata', metadata, 'Compression', 'LZW' );
    answer = true;
end

