function preprocessing(imgs_paths,flag_vector,aspect,min_obj_size)
	
	current_path = which(mfilename);
	[current_path, filename, extension] = fileparts( current_path );
	cd(current_path);

	%loop each image
	factor_matrices={};
	mkdir temp
	for f=1:length(imgs_paths)
		img_path=imgs_paths{f};
		reader = bfGetReader(img_path);
		omeMeta = reader.getMetadataStore();
		z_size = omeMeta.getPixelsSizeZ(0).getValue();

		%loop each channel
		puncta_or_threshold_flags=[];
		other_channel_imgs={};
		dist_trans_matrices={};
		for c=1:length(flag_vector)
			flag=flag_vector{c};
			img=imreadBF(img_path,1:z_size,1,c);
			
			if ~isempty(flag)
				if strcmp(flag{2},'model')
					if strcmp(flag{1},'puncta')
						[~,centers]=find_punta(img,'puncta');
					end
					if strcmp(flag{1},'objects')
						[~,centers]=find_objects(img,'objects_center',min_obj_size);
					end
				end
				if strcmp(flag{2},'factor')
					if strcmp(flag{1},'puncta')
						[puncta_img,~]=find_punta(img,'puncta_img');
						disp('transformation matrix caculation start...');
						dist_trans_matrices{end+1}=bwdistsc(puncta_img,aspect);
						save3Dmatrix(dist_trans_matrices{end},img_path,length(dist_trans_matrices));
						disp('transformation matrix caculation complete');
					end
					if strcmp(flag{1},'objects')
						[objects_img,~]=find_objects(img,'objects_img',min_obj_size);
						disp('transformation matrix caculation start...');
						dist_trans_matrices{end+1}=bwdistsc(objects_img,aspect);
						save3Dmatrix(dist_trans_matrices{end},img_path,length(dist_trans_matrices));
						disp('transformation matrix caculation complete');
					end
				end
			else
				continue;
			end
		end

		%save centers into txt file
		centers_num=size(centers,2);
		centers2=zeros(centers_num,3);
		for i=1:centers_num
			centers2(i,:)=centers{i};
		end
		pre_name=strsplit(img_path,filesep);
		centers_filename = ['./temp/' pre_name{end} '__centers.txt'];
    	save(centers_filename,'centers2','-ASCII','-double');
		
		%factor_matrix
		dist_trans_matrices_num=size(dist_trans_matrices,2);
		factor_matrix=zeros(centers_num,dist_trans_matrices_num);
	    for i=1:centers_num
	        for j=1:dist_trans_matrices_num
	            factor_matrix(i,j)=dist_trans_matrices{j}(centers{i}(1),centers{i}(2),centers{i}(3));
	        end
	    end

	    %save factor_matrix into txt file
		factor_matrix_filename = ['./temp/' pre_name{end} '__factor_matrix.txt'];
    	save(factor_matrix_filename,'factor_matrix','-ASCII','-double');
	end
end

function save3Dmatrix(m,img_path,dm_n)
	% sl3-AR3-AR4-NeuN-MBP-VDAC-40x004.nd2__dist_trans_matrix_2_z9.txt
	s=strsplit(img_path,filesep);
	
	for z=1:size(m,3)
		filename=['./temp/' s{end} '__dist_trans_matrix_' int2str(dm_n) '_z' int2str(z) '.txt'];
		disp(filename);
		A=m(:,:,z);
		save(filename,'A','-ASCII','-double');
	end
end