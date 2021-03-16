function tiff2ometiff_2D()
    cd('~/CellOrganizor/build/cellorganizer3/');
    setup;
    pattern = {'LAM' 'Mit' 'Nuc' 'TfR'};
    for iii=1:4
        image_directory = ['./images/HeLa/2D/' pattern{iii} filesep 'crop'];
        system(['mkdir ./ometiff_image/hela/hela_2D_' lower(pattern{iii})]);
        number_of_files = length( dir([image_directory filesep 'cell*.tif']) );
        for i=1:number_of_files
            disp(['Creating 2D image ' num2str(i)]);
            image1_filename=['./images/HeLa/2D/' pattern{iii} filesep 'orgdna' filesep 'cell' num2str(i) '.tif'];
            image2_filename=['./images/HeLa/2D/' pattern{iii} filesep 'orgcell' filesep 'cell' num2str(i) '.tif'];
            image3_filename=['./images/HeLa/2D/' pattern{iii} filesep 'orgprot' filesep 'cell' num2str(i) '.tif'];
            image4_filename=['./images/HeLa/2D/' pattern{iii} filesep 'crop_uint8' filesep 'cell' num2str(i) '.tif'];
            %%check image channel consistence
            image1=tif2img(image1_filename);
            image2=tif2img(image2_filename);
            image3=tif2img(image3_filename);
            image4=tif2img(image4_filename);
            v=[size(image1,3),size(image2,3),size(image3,3),size(image4,3)];
            if ~all(v == v(1))
                error('image channel not consistent.');
            end
            %save to ometiff
            images = {image1_filename,image2_filename,image3_filename,image4_filename};
            pre=['hela_2D_' lower(pattern{iii}) '_cell_' num2str(i)];
            filename = ['./ometiff_image/hela/hela_2D_' lower(pattern{iii}) filesep  pre '.ome.tif'];
            parameters.PhysicalSizeX = 0.0490;
            parameters.PhysicalSizeY = 0.0490;
            parameters.PhysicalSizeZ = 0.2000;
            parameters.list_of_channel_labels = {[pre '_nucleus'], [pre '_cell'], [pre '_protein'], [pre '_crop']};
            disp( ['Saving image as ' filename ]);
            tif2ometiff( images, filename, parameters  );
        end
    end
end