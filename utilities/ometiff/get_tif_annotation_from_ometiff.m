function get_tif_annotation_from_ometiff( directory )
% GET_TIF_ANNOTATION_FROM_OMETIFF Retrieve tif header from OME.TIFF

% Xin Lu
%
% Copyright (C) 2018-2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

image_paths=ml_ls([directory filesep '*.ome.tif']);
for i = 1:length(image_paths)
    tiff=Tiff(image_paths{i});
    ImageDescription=tiff.getTag('ImageDescription');
    newStr = extractBetween(ImageDescription,'<UUID FileName="','">'); %extract original ometiff name
    original_ometiff_name=newStr{1};
    original_ometiff_name_=strrep(original_ometiff_name,'.ome.tif','');
    
    if ~exist( [ pwd filesep original_ometiff_name] )
        if ~strcmp(image_paths{i},original_ometiff_name)
            system(['cp "' image_paths{i} '" "' pwd filesep original_ometiff_name '"']);
        end
        % image_path=[ometiff_img_path '/' char(img_file_name_pre(i)) '.ome.tif'];
        % GFP_foldername=['./temp_for_ometiff/"' strrep(char(img_file_name_pre(i)),'_','/') '"/GFP'];
        GFP_foldername=['./temp_for_ometiff/"' strrep(char(original_ometiff_name_),'_','/') '"/GFP'];
        image_path=original_ometiff_name;
        
        reader = bfGetReader(image_path);
        omeMeta = reader.getMetadataStore();
        
        %extract tif from ometiff
        system(['mkdir -p ' GFP_foldername]);
        x_size = omeMeta.getPixelsSizeX(0).getValue();
        y_size = omeMeta.getPixelsSizeY(0).getValue();
        z_size = omeMeta.getPixelsSizeZ(0).getValue();
        t_size = omeMeta.getPixelsSizeT(0).getValue();
        whole_img = tif2img(image_path);
        img_DIC=zeros(y_size,x_size,t_size);
        for t = 1:t_size
            % t=40
            fprintf('Extracting timepoint %d in %s ...\n',t,char(original_ometiff_name_));
            img = zeros(y_size,x_size,z_size-1);
            for z = 1:z_size-1
                img(:,:,z) = whole_img(:,:,z*t_size-t_size+t);
            end
            tt=int2str(t);
            if ismember(t,1:9)
                tt=['0' tt];
            end
            cd([pwd filesep strrep(GFP_foldername(2:end),'"','')]);
            img2tif(img,['GFP (Timepoint ' tt ').tif']);
            cd('../../../../');
            img_DIC(:,:,t)=whole_img(:,:,t_size*(z_size-1)+t);
        end
        cd([pwd strrep(GFP_foldername(2:end-4),'"','')]);
        img2tif(img_DIC,'DIC.tif');
        cd('../../../');
        
        %extract annotation from ometiff
        system(['mkdir -p ./temp_for_ometiff/annotations']);
        map_annotation=char(omeMeta.getMapAnnotationValue(0));
        map_split=strsplit(map_annotation,', MapPair<annotation');
        if strcmp(original_ometiff_name_,'5C.C7 LAT 15 03 02_run 1')
            annot_file_name_pre='coordinates 5C.C7 LAT 02 03 15.csv';
        else
            str1=strsplit(original_ometiff_name_,'_');
            str2=strsplit(str1{1},' ');
            % '5C.C7 LAT 15 03 16_run 1'
            % 'coordinates 5C.C7 LAT 16 03 15 run 1.csv'
            annot_file_name_pre=[ 'coordinates 5C.C7 LAT ' str2{5} ' ' str2{4} ' 15 ' str1{2} '.csv'];
        end
        csv_filename=['./temp_for_ometiff/annotations/' char(annot_file_name_pre)];
        if exist(csv_filename,'file') == 2
            system(['rm "' csv_filename '"']);
        end
        fileID = fopen(csv_filename,'a');
        for j=1:length(map_split)
            m2=char(map_split(j));
            m3=strsplit(m2,'\n');
            m4=char(m3(2));
            m5=strtrim(m4);
            fprintf(fileID,'%s\n',m5);
        end
        fclose(fileID);
    else
       disp(['Image ' original_ometiff_name ' exists. Skipping file']); 
    end
end
end