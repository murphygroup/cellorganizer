function tiff2ometiff_3D()
    cd('~/CellOrganizor/build/cellorganizer3/');
    setup;
    pattern = {'LAM' 'Mit' 'Nuc' 'TfR'};
    tp={'raw' 'processed'};
    for iii=1:4
        for ii=1:2
            image_directory = ['./images/HeLa/3D/' tp{ii}];
            system(['mkdir ./ometiff_image/hela/hela_3D_' tp{ii} '_' lower(pattern{iii})]);
            number_of_files = length( dir([image_directory filesep pattern{iii} '_cell*_ch0_t1.tif']) );
            for i=1:number_of_files
                disp(['Creating 3D image ' num2str(i)]);
                %process mask
                img_0 = tif2img( ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{ii} filesep pattern{iii} '_cell' num2str(i) '_ch0_t1.tif'] );
                img_mask = tif2img( ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{ii} filesep pattern{iii} '_cell' num2str(i) '_mask_t1.tif'] );
                if size(img_mask,3)>1
                    img_mask=img_mask(:,:,1)+img_mask(:,:,2);
                end
                img_mask = repmat(img_mask,1,1,size(img_0,3));
                img2tif(img_mask, ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{ii} filesep pattern{iii} '_cell' num2str(i) '_mask_update_t1.tif']);

                images = {[image_directory filesep pattern{iii} '_cell' num2str(i) '_ch0_t1.tif'], ...
                    [image_directory filesep pattern{iii} '_cell' num2str(i) '_ch1_t1.tif'], ...
                    [image_directory filesep pattern{iii} '_cell' num2str(i) '_ch2_t1.tif'], ...
                    [image_directory filesep pattern{iii} '_cell' num2str(i) '_mask_update_t1.tif']};
                pre=['hela_3D_' tp{ii} '_' lower(pattern{iii}) '_cell_' num2str(i)];
                filename = ['./ometiff_image/hela/hela_3D_' tp{ii} '_' lower(pattern{iii}) filesep  pre '.ome.tif'];
                parameters.PhysicalSizeX = 0.0490;
                parameters.PhysicalSizeY = 0.0490;
                parameters.PhysicalSizeZ = 0.2000;
                parameters.list_of_channel_labels = {[pre '_nucleus'], [pre '_cell'], [pre '_protein'], [pre '_mask']};
                disp( ['Saving image as ' filename ]);
                tif2ometiff( images, filename, parameters  );
            end
        end
    end

    % i=34;
    % img_0 = tif2img( ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{1} filesep pattern{1} '_cell' num2str(i) '_ch0_t1.tif'] );
    % img_mask = tif2img( ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{1} filesep pattern{1} '_cell' num2str(i) '_mask_t1.tif'] );
    % if size(img_mask,3)>1
    %     img_mask=img_mask(:,:,1)+img_mask(:,:,2);
    %     111
    % end
    % img_mask = repmat(img_mask,1,1,size(img_0,3));
    % img2tif(img_mask, ['~/CellOrganizor/build/cellorganizer3/images/HeLa/3D/' tp{1} filesep pattern{1} '_cell' num2str(i) '_mask_update_t1.tif']);
end