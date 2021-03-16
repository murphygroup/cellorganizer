function get_FactorMatrix(imgs_paths,masks_paths,intermediate_output_path,flag_vector,mc,options)
	others_output_path=[intermediate_output_path filesep 'others'];
	if ~exist(others_output_path)
		system(['mkdir -p ' others_output_path]);
	end
	factor_channels=get_factor_channels(flag_vector);
	%loop each image
	for f=1:length(imgs_paths)
		img_path=imgs_paths{f};
		if ~is_ome_tiff(img_path)
			z_size=get_z_size(img_path);
		end
		pre_name=strsplit(img_path,filesep);
		%for mc
		if is_ome_tiff(img_path)
			img=OME_loadchannel(img_path,mc);
		else
			img=imreadBF(img_path,1:z_size,1,mc);
		end
		if ~isempty(masks_paths)
			[labeled_mask,maskIDs]=get_labeled_mask(masks_paths{f},options.mask_inverted_color_flag);
			[masked_imgs,ROIIDs]=get_masked_imgs(img,labeled_mask,maskIDs);
        else
        	masked_imgs={img};
        	ROIIDs={0};
		end
		%loop each masked_imgs
		for m=1:length(masked_imgs)
            centers_filename=[others_output_path filesep 'centers__' pre_name{end} '-' int2str(ROIIDs{m}) '.txt']; %haven't provided mc info since no need to
			[~,centers]=find_punta(img,'puncta');
			save_centers(centers,centers_filename);
			dist_trans_matrices=read_dist_trans_matrices(factor_channels,intermediate_output_path,pre_name,ROIIDs{m});
			centers=dlmread(centers_filename); %read back centers
			factor_matrix=make_factor_matrix(dist_trans_matrices,centers);
			save_factor_matrix(others_output_path,pre_name,ROIIDs{m},factor_matrix);
			if options.debug
				for_debug(intermediate_output_path,img_path,factor_matrix,pre_name,mc,ROIIDs{m},centers);
			end
		end
	end
end

function for_debug(intermediate_output_path,img_path,factor_matrix,pre_name,mc,ROIID,centers)
	debug_output_path=[intermediate_output_path filesep 'debug'];
	if ~exist(debug_output_path)
		system(['mkdir -p ' debug_output_path]);
	end
	[x_size,y_size]=get_size(img_path);
	for j=1:size(factor_matrix,2)
		filename=[debug_output_path filesep 'factor_matrix__' pre_name{end} '-' int2str(ROIID) '-' int2str(mc) '-' int2str(j) '_debug.jpg'];
		img_debug=zeros(x_size,y_size);
		for jj=1:size(factor_matrix,1)
			img_debug(centers(jj,1),centers(jj,2))=factor_matrix(jj,j);
		end
		imwrite(img_debug,filename);
	end
end

function [x_size,y_size]=get_size(img_path)
	reader = bfGetReader(img_path);
	omeMeta = reader.getMetadataStore();
	x_size = omeMeta.getPixelsSizeX(0).getValue();
	y_size = omeMeta.getPixelsSizeY(0).getValue();
end

function save_centers(centers,centers_filename)
	centers_num=size(centers,2);
	centers2=zeros(centers_num,3);
	for i=1:centers_num
		centers2(i,:)=centers{i};
	end
	disp(['saving centers...']);
	save(centers_filename,'centers2','-ASCII','-double');
end

function dist_trans_matrices=read_dist_trans_matrices(factor_channels,intermediate_output_path,pre_name,ROIID)
	disp(['reading back dist_trans_matrices...']);
	dist_trans_matrices={};
	for c_ind=1:length(factor_channels)
		data=load([intermediate_output_path filesep 'dist_trans_matrix' filesep 'dist_trans_matrix__' pre_name{end} '-' int2str(ROIID) '-' int2str(factor_channels(c_ind)) '.mat']);
		dist_trans_matrices{end+1}=data.m;
	end
end

function factor_matrix=make_factor_matrix(dist_trans_matrices,centers)
	dist_trans_matrices_num=size(dist_trans_matrices,2);
	centers_num=size(centers,1);
	factor_matrix=zeros(centers_num,dist_trans_matrices_num);
    for i=1:centers_num
        for j=1:dist_trans_matrices_num
            % factor_matrix(i,j)=dist_trans_matrices{j}(centers(i,1),centers(i,2),centers(i,3));
            factor_matrix(i,j)=dist_trans_matrices{j}(centers(i,2),centers(i,1),centers(i,3));
        end
    end
end

function save_factor_matrix(others_output_path,pre_name,ROIID,factor_matrix)
	factor_matrix_filename = [others_output_path filesep 'factor_matrix__' pre_name{end} '-' int2str(ROIID) '.txt'];
	disp(['saving factor matrix...']);
	save(factor_matrix_filename,'factor_matrix','-ASCII','-double');
end

function factor_channels=get_factor_channels(flag_vector);
	factor_channels=[];
	for i=1:length(flag_vector)
		flag=flag_vector{i};
		if ~isempty(flag) && any(strcmp(flag{2},'factor'))
			factor_channels(end+1)=i;
		end
	end
end