function [big_matrix, name_matrix] = ...
    tz_process3_ds( dirname1, dirname2, dirname3, task, ...
    ext, threshold, medfilter, level)
%TZ_PROCESS3_DS Process and analyze drug images.
%   TZ_PROCESS3_DS(DIRNAME1,DIRNAME2,DIRNAME3,TASK,EXT,THRESHOLD,MEDFILTER)
%   does one of the three tasks, projecting, making mask and calculating 
%   features by specifying TASK.
%-------------------------------------------------------
%  TASK          DIRNAME1      DIRNAME2        DIRNAME3
%'make_projs'   drugimages    projimages         empty
%'make_masks'   projimages    maskimages         empty
%'calc_feats'   drugimages    maskimages       features
%-------------------------------------------------------
%
%   EXT is for extension of the image files. Any pixel value greater
%   than THRESHOLD will be set to 0. If MEDFILTER is 1, all images 
%   will be preprocessed by a median filter. If LEVEL is 'cell', it is
%   for cell level. If it is 'object', it is for oject level processing.
%   
%   [BIG_MATRIX,NAME_MATRIX] = TZ_PROCESS3_DS(DIRNAME1,DIRNAME2,DIRNAME3,
%   'calc_feats',...) also returns features and names for each feature
%   matrix.

%   17-Sep-2005 Initial write T. Zhao
%   ??-???-???? Modified from Meel's function T. Zhao
%   06-MAR-2003 Modified T. Zhao
%   26-MAR-2003 Modified T. Zhao
%   01-APR-2003 Modified T. Zhao
%   18-MAY-2003 Modified T. Zhao
%   15-MAR-2005 Modified T. Zhao
%       - change xc_loadimage to tz_loadimage
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 7
    error('7 or 8 arguments are required')
end

if ~exist('level','var')
    level='cell';
end

prot_list = ml_dir(dirname1)
prot_list = prot_list(3:end)
big_matrix = {};
for prot_no = 1:length(prot_list)
    prot_name = prot_list{prot_no}
    prot_fullpath = [dirname1 '/' prot_name];
    prot_fullpath2 = [dirname2 '/' prot_name];
    prot_fullpath3 = [dirname3 '/' prot_name];
        
    switch( task)
    case {'make_projs','make_projs_sum','make_masks'},
        prot_fullpath2 = [dirname2 '/' prot_name];
        if ~exist(dirname2,'dir')
            make_dir(dirname2);
        end
        make_dir( prot_fullpath2);
    case {'calc_feats','calc_feats_texture','calc_newfeats'},
        if ~exist(dirname3,'dir')
            make_dir(dirname3);
        end
        make_dir( prot_fullpath3);
    otherwise
        error('Incorrect task name!');
    end
    drug_list = ml_dir(prot_fullpath);
    drug_list = drug_list(3:end);
    ncondit = length(drug_list);
    ndrugs = ncondit/2;
    for drug_no = 1:ncondit
        drug_name = drug_list{drug_no};
        drug_fullpath = [prot_fullpath '/' drug_name];
        drug_fullpath2 = [prot_fullpath2 '/' drug_name];
        drug_fullpath3 = [prot_fullpath3 '/' drug_name];    
        switch( task)
        case {'make_projs','make_projs_sum','make_masks'},
            make_dir( drug_fullpath2);
        case {'calc_feats','calc_feats_texture','calc_newfeats'},
            make_dir( drug_fullpath3);
        otherwise
            error('Incorrect task name!');
        end
        features = process_cond( drug_fullpath, drug_fullpath2, drug_fullpath3, task, ext, threshold, medfilter, level);
        % assumes that all the with drug conditions come first and then all of the 
        % no drug conditions
        if (drug_no > ndrugs)
            big_matrix{1, drug_no-ndrugs, prot_no} = features;
            name_matrix{1, drug_no-ndrugs, prot_no} = {prot_name, drug_name};
        else
            big_matrix{2, drug_no, prot_no} = features;
            name_matrix{2, drug_no, prot_no} = {prot_name, drug_name};
        end
        %save('big_matrix_3.mat','big_matrix');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function feature_matrix = process_cond( dirname1, dirname2, dirname3, task, ext, threshold, medfilter, level)

dirname1
image_names = ml_dir(dirname1);
% image_names
image_names = tz_cleandirs(image_names(3:end))

no_of_images = length(image_names);
feature_matrix = [];
for image_no = 1:no_of_images
    image_name = image_names{image_no}
    
    %Do not process time series images
    %Matlab 5.3 does not support continue
    process=1;
    
    if length(image_name)>=7
        if all(image_name(1:7)=='timeser')
            process=0
        end
    end
    
%    process
%    pause

    if process==1
        image_fullname = [dirname1 '/' image_name]
        switch( task)
        case 'calc_newfeats'
            feat_fullpath = [dirname3 '/' image_name];
            if( ~exist( [feat_fullpath '.mat']) )
                % Read image
                image = ml_loadimage( image_fullname, ext, threshold);
                
                no_of_slices = size( image, 3);
                if medfilter == 1
                    %denoise
                    'medfilter'
                    for slice_no = 1:no_of_slices
                        slice = image(:,:,slice_no);
                        %denoise
                        image(:,:,slice_no)=medfilt2(slice,[3,3]);
                    end
                end
                
                % Mask the images
                mask_img_fullpath = [dirname2 '/' image_name '.mat']
                load ( mask_img_fullpath);
                
                imagesize=size(image)
                
                % Compute features
                switch( level)
                case 'cell'
                    scale=[0.1075 0.5];
                    tratio = repmat(0.5,[1,length(scale)])./scale;
                    tgray=[];
                    [names, image_features, slfnames] = ml_3dfeatset( image, 'SLF19', ...
						  mask,[],tratio,tgray,[1 1 203/48.8],'yesbgsub','rc');
                    
                case 'object'
                    error('Sorry, object level is unavailable');
                otherwise
                    error('Incorrect level name!');
                end
                % Save features for this image
                save( feat_fullpath,'image_features');
            else
                load( feat_fullpath);
            end
            % Assemble feature matrix for all images
            
            feature_matrix = [ feature_matrix; image_features];
        case 'calc_feats',
            feat_fullpath = [dirname3 '/' image_name]
            if( ~exist( [feat_fullpath '.mat']) )
                % Read image
                image = ml_loadimage( image_fullname, ext, threshold);
                
                no_of_slices = size( image, 3);
                if medfilter == 1
                    %denoise
                    'medfilter'
                    for slice_no = 1:no_of_slices
                        slice = image(:,:,slice_no);
                        %denoise
                        image(:,:,slice_no)=medfilt2(slice,[3,3]);
                    end
                end
                
                % Convert image to double and substract the background
                image=double(image);
                max_pixel = max(image(:));
                image=image/(max_pixel/256);
                image=uint8(image);
                image = ml_3dbgsub(image);
                % Mask the images
                mask_img_fullpath = [dirname2 '/' image_name '.mat']
                load ( mask_img_fullpath);
                
                for slice_no = 1:no_of_slices
                    slice = image(:,:,slice_no);
                    slice(find(mask==0))=0;
                    image(:,:,slice_no) = slice;
                end
                imagesize=size(image)

                % Compute features
                switch( level)
                case 'cell'
                    image_features = ml_3dfeatures( image, [], [], [0.1075 0.1075 0.5]);
                case 'object'
                    image_features = tz_calobj( image, [], [], [0.1075 0.1075 0.5]);
                otherwise
                    error('Incorrect level name!');
                end
                % Save features for this image
                save( feat_fullpath,'image_features');
            else
                load( feat_fullpath);
            end
            % Assemble feature matrix for all images
            
            feature_matrix = [ feature_matrix; image_features];
            
        case 'calc_feats_texture',
            feat_fullpath = [dirname3 '/' image_name]
            if( ~exist( [feat_fullpath '.mat']) )
                             
               
          
                mask_img_fullpath = [dirname2 '/' image_name '.mat']
                %load ( mask_img_fullpath);
                image_fullname
                mask_img_fullpath
                image_features = tz_texture_3d_ds( image_fullname, mask_img_fullpath,ext, threshold,medfilter);

               
                % Save features for this image
                save( feat_fullpath,'image_features');
            else
                load( feat_fullpath);
            end
            % Assemble feature matrix for all images
            
            feature_matrix = [ feature_matrix; image_features];
        case 'make_projs',
            proj_img_fullpath = [dirname2 '/' image_name '.mat']
            if( ~exist( proj_img_fullpath))
                % Read image
                image = ml_loadimage( image_fullname, ext, threshold); 
                if isempty(image)
                    return;
                end
                % Test for bad thresholding
                %nr_thresh_pixels = length(find(image == threshold));
                %f = fopen( '/home/sl/thresh_values.txt', 'a');
                %fprintf( f, [image_fullname num2str( nr_thresh_pixels) '\n']);
                %fclose(f);
                % Create projection image
                proj_image = max( image, [], 3);
                save( proj_img_fullpath, 'proj_image');
            end
            
        case 'make_projs_sum',
            proj_img_fullpath = [dirname2 '/' image_name '.mat']
            if( ~exist( proj_img_fullpath))
                % Read image
                image = ml_loadimage( image_fullname, ext, threshold); 
                
                % Test for bad thresholding
                %nr_thresh_pixels = length(find(image == threshold));
                %f = fopen( '/home/sl/thresh_values.txt', 'a');
                %fprintf( f, [image_fullname num2str( nr_thresh_pixels) '\n']);
                %fclose(f);
                % Create projection image
                proj_image = sum( image, 3);
                save( proj_img_fullpath, 'proj_image');
            end           
            
        case 'make_masks',
            mask_img_fullpath = [dirname2 '/' image_name];
            if( ~exist( mask_img_fullpath))
                % Read projection image
                load( image_fullname);
                image_size = size( proj_image);
                imagesc( proj_image);
                truesize( image_size/1.5);
                title(mask_img_fullpath);
                mask = roipoly;
                save( mask_img_fullpath, 'mask') ;        
            end
            
        otherwise,
            error('Incorrect task name!');
        end
        
    end  % end of " if image_name(1:7) ~= 'timeser' "
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function make_dir( dirname)

if( ~exist( dirname, 'dir'))
    command = ['mkdir ' dirname]
    status = unix( command);
    if( status ~= 0)
        error(['Cannot create directory: ' dirname]);
    end
end

