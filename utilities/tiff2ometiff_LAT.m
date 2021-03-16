function tiff2ometiff_LAT()
    cd('~/CellOrganizor/build/cellorganizer3/');
    setup;
    csv_dir='./images/LAT/annotations';
    csv=dir([csv_dir filesep '*.csv']);
    for i=1:length(csv);
        % i=1;%%%%%
        disp(['csv: ' csv(i).name]);
        if ~isempty(findstr(csv(i).name,' '))
            nn1=strrep(csv(i).name,' ','_');
            nn2=strrep(csv(i).name,' ','\ ');
        end
        system(['mkdir -p ./ometiff_image/lat/' nn1]);

        [stat,cmd]=system(['grep ").tif" ' csv_dir filesep nn2]);
        cmd=strsplit(cmd,'\n');

        for ii=1:length(cmd)-1
            % ii=1;%%%%%
            file=strsplit(cmd{ii},',');
            full_name=file{1};
            full_name=strrep(full_name,'..','./images/LAT');
            images{1}=full_name;
            parameters.list_of_channel_labels{1}=file{1};
            parameters.PhysicalSizeX = 0.0490;
            parameters.PhysicalSizeY = 0.0490;
            parameters.PhysicalSizeZ = 0.2000;
            fname=strsplit(full_name,'/');
            tifname=fname{end};
            tifname=tifname(1:end-4);
            if ~isempty(findstr(tifname,' '))
                tifname=strrep(tifname,' ','_');
            end
            tif2ometiff(images,['./ometiff_image/lat/' nn1 filesep tifname '.ome.tif'],parameters);
        end
    end
end