function answer = get_all_distance_transformation_matrix(imgs_paths,masks_paths,intermediate_output_path,flag_vector,options)
% get_all_distance_transformation_matrix generates all distance transformation matrices for all factor channels of all images. The resulting distance transformation matrices will be saved in intermediate_output_path. For example, if the flag_vector argument is set as {{},{'object',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}}}, for each image, the function will generate three distance transformation matrices for channel 2, 3, 4.
%
% Input 
% -----
% imgs_paths: A cell of strings that constains images' paths.
% masks_paths: A cell of strings that constains masks' paths. The order in masks_paths should match the order in imgs_paths. It can be set as empty.
% intermediate_output_path: The path that saves the intermediate files, which will be used as the path to save the resulting distance transformation matrices for this function.
% flag_vector: A three dimension cell containing flags. It can contains
%				0 or 1 from
%					puncta = find puncta for this channel
%					objects = find objects for this channel
%                   objects_local = find small objects for this channel (for mitochondria)
%					multiple_objects = find multiple objects for this channel
%					intensity = remain image intensity for this channel
%				0, 1 or 2 from
%					model = build model for the centers after preprocessing this channel
%					factor = use as factor after preprocessing and distance transform
%				e.g., {{},{'object',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}}} means
%				1. build a model of puncta 1 (channel 3) using one factor: distance to object found using channel 2; 
%				2. build a model of puncta 1 (channel 3) using one factor: distance to puncta found using channel 4; 
%				3. build a model of puncta 1 (channel 3) using two factors: distance to object found using channel 2 and distance to puncta found using channel 4; 
%				4. build a model of puncta 2 (channel 4) using one factor: distance to object found using channel 2; 
%				5. build a model of puncta 2 (channel 4) using one factor: distance to puncta found using channel 3; 
%				6. build a model of puncta 2 (channel 4) using two factors: distance to object found using channel 2 and distance to puncta found using channel 3. 
%				Note: only puncta channel can be set as model channel. Channels that have empty flags mean they will be ignored. This argument must be the same flag_vector as in function get_all_models. In this function, flag_vector will be used only for generating all distance transformation matrices for all factor channels of all images. For the given example, for each image, the function will generate three distance transformation matrices for channel 2, 3, 4.
% options.aspect: Aspect ratio for each image dimension.
% options.min_obj_size: The threshold of size that is used for determining an object when finding an object.
% options.debug: A boolean value indicates whether visualize the distance transformation matrix. The resulting figures will save in intermediate_output_path/debug folder.
% options.mask_inverted_color_flag: A boolean value indicates whether mask color need to invert.
%
% Output
% ------
% answer: A boolean indicates whether the function successfully finished or not.

% Author: Xin Lu
%
% Copyright (C) 2012-2019 Murphy Lab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
answer=1;
% try

	flag_vector=get_flag_vector_slim(flag_vector);
	dist_trans_matrix_output_path=[intermediate_output_path filesep 'dist_trans_matrix'];
	if ~exist(dist_trans_matrix_output_path)
		system(['mkdir -p ' dist_trans_matrix_output_path]);
	end
	%loop each image
	for f=1:length(imgs_paths)
	% f=1;
		img_path=imgs_paths{f};
		if ~is_ome_tiff(img_path)
			z_size=get_z_size(img_path);
		end
		if ~isempty(masks_paths)
			[labeled_mask,maskIDs]=get_labeled_mask(masks_paths{f},options.mask_inverted_color_flag);
			% save_cell_boundary(labeled_mask,maskIDs);
		end
		%loop each channel
		for c=1:length(flag_vector)
			flag=flag_vector{c};
			% if c==2
			% 	2
			% end
			% if c==1
			% 	continue
			% end
			if ~isempty(flag)
				if is_ome_tiff(img_path)
					img=OME_loadchannel(img_path,c);
				else
					img=imreadBF(img_path,1:z_size,1,c);
				end
				if ~isempty(masks_paths)
					[masked_imgs,ROIIDs]=get_masked_imgs(img,labeled_mask,maskIDs);
		        else
		        	masked_imgs={img};
		        	ROIIDs={0};
				end
				%loop each masked_imgs
				for m=1:length(masked_imgs)
					if strcmp(flag{1},'puncta')
						[puncta_img,~]=find_punta(masked_imgs{m},'puncta_img');
						get_distance_transformation_matrix(puncta_img,img_path,intermediate_output_path,options.aspect,c,ROIIDs{m},options);
					end
					if strcmp(flag{1},'objects')
						[objects_img,~]=find_objects(masked_imgs{m},'objects_img',options.min_obj_size);
						get_distance_transformation_matrix(objects_img,img_path,intermediate_output_path,options.aspect,c,ROIIDs{m},options);
                    end
                    if strcmp(flag{1},'objects_local') % options for mitochondria
						[objects_img,~]=find_objects_local(masked_imgs{m},'objects_img',options.min_obj_size);
						get_distance_transformation_matrix(centers,img_path,intermediate_output_path,options.aspect,c,ROIIDs{m},options);
					end
					if strcmp(flag{1},'boundaries')
						[boundaries_img]=find_boundaries(masked_imgs{m});
						get_distance_transformation_matrix(boundaries_img,img_path,intermediate_output_path,options.aspect,c,ROIIDs{m},options);
					end
					if strcmp(flag{1},'nucleus')
						[nucleus_img,~]=find_nucleus(masked_imgs{m},'nucleus_img',options.min_obj_size);
						get_distance_transformation_matrix(nucleus_img,img_path,intermediate_output_path,options.aspect,c,ROIIDs{m},options);
					end
				end
			end
			% break
		end
		% break
	end

% catch
% 	answer=0;
% 	warning('get_all_distance_transformation_matrix has error');
% end
end

function flag_vector_slim=get_flag_vector_slim(flag_vector)
flag_vector_slim={};
	for i=1:length(flag_vector)
		flag=flag_vector{i};
		if isempty(flag)
			flag_vector_slim{end+1}={};
		elseif any(strcmp(flag{2},'factor'))
			flag_vector_slim{end+1}={flag{1},'factor'};
		else
			flag_vector_slim{end+1}={};
		end
	end
end

function get_distance_transformation_matrix(processed_img,img_path,intermediate_output_path,aspect,c,ROIID,options)
	disp([img_path '--channel:' num2str(c) '--ROIID:' num2str(ROIID) ' transformation matrix caculation start...']);
	dist_trans_matrix=bwdistsc(processed_img,aspect);
	save3Dmatrix(dist_trans_matrix,img_path,c,ROIID,intermediate_output_path);
	if options.debug
		debug_output_path=[intermediate_output_path filesep 'debug'];
		if ~exist(debug_output_path)
			system(['mkdir -p ' debug_output_path]);
		end
		s=strsplit(img_path,filesep);
		filename=[debug_output_path filesep 'dist_trans_matrix__' s{end} '-' int2str(ROIID) '-' int2str(c) '_debug.jpg'];
		imwrite(z_projection(dist_trans_matrix),filename);
	end
end

function save3Dmatrix(m,img_path,dm_n,ROIID,intermediate_output_path)
	s=strsplit(img_path,filesep);
	filename=[intermediate_output_path filesep 'dist_trans_matrix' filesep 'dist_trans_matrix__' s{end} '-' int2str(ROIID) '-' int2str(dm_n) '.mat'];
	disp(['saving ' filename]);
	save(filename,'m');
end

function projection_img=z_projection(img)
	projection_img=zeros(size(img,1),size(img,2));
	z_size=size(img,3);
	for z=1:z_size
		projection_img=projection_img+img(:,:,z);
	end
	projection_img=projection_img/z_size;
end