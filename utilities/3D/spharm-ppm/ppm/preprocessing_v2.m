function preprocessing(imgs_paths,flag_vector,aspect,min_obj_size,model_flag)
	
	current_path = which(mfilename);
	[current_path, filename, extension] = fileparts( current_path );
	cd(current_path);
	
	%loop each image
	factor_matrices=[];
	if ~exist( './temp/dist_trans_matrix' )
		mkdir ./temp/dist_trans_matrix
	end
	if ~exist( './temp/others' )
		mkdir ./temp/others
	end
	for f=1:length(imgs_paths)
	% f=1;
		img_path=imgs_paths{f};
		puncta_or_threshold_flags=[];
		other_channel_imgs={};
		dist_trans_matrices={};
		factor_flag=[];
		%loop each channel
		for c=1:length(flag_vector)
			% c
			flag=flag_vector{c};
			if ~model_flag
				if c==1
					reader = bfGetReader(img_path);
					omeMeta = reader.getMetadataStore();
					z_size = omeMeta.getPixelsSizeZ(0).getValue();
				end
				img=imreadBF(img_path,1:z_size,1,c);
				if strcmp(flag{1},'puncta')
					[puncta_img,~]=find_punta(img,'puncta_img');
					disp('transformation matrix caculation start...');
					dist_trans_matrices{end+1}=bwdistsc(puncta_img,aspect);
					save3Dmatrix(dist_trans_matrices{end},img_path,c);
				end
				if strcmp(flag{1},'objects')
					[objects_img,~]=find_objects(img,'objects_img',min_obj_size);
					disp('transformation matrix caculation start...');
					dist_trans_matrices{end+1}=bwdistsc(objects_img,aspect);
					save3Dmatrix(dist_trans_matrices{end},img_path,c);
				end
			else
				if ~isempty(flag)
					if strcmp(flag{2},'model')
						% celldisp(flag);
						factor_flag=[factor_flag,0];
						if model_flag
							pre_name=strsplit(img_path,filesep);
							centers_filename = ['./temp/others/' pre_name{end} '__centers.txt'];
							if exist(centers_filename)~=2
								reader = bfGetReader(img_path);
								omeMeta = reader.getMetadataStore();
								z_size = omeMeta.getPixelsSizeZ(0).getValue();
								img=imreadBF(img_path,1:z_size,1,c);
								if strcmp(flag{1},'puncta')
									[~,centers]=find_punta(img,'puncta');
								end
								if strcmp(flag{1},'objects')
									[~,centers]=find_objects(img,'objects_center',min_obj_size);
								end
								%save centers into txt file
								centers_num=size(centers,2);
								centers2=zeros(centers_num,3);
								for i=1:centers_num
									centers2(i,:)=centers{i};
								end
								disp(['saving centers...']);
						    	save(centers_filename,'centers2','-ASCII','-double');
						    end
						end
					end
					if strcmp(flag{2},'factor')
						factor_flag=[factor_flag,1];
					end
				else
					factor_flag=[factor_flag,0];
				end
			end
		end

		if model_flag
			%read back dist_trans_matrices
			disp(['reading back dist_trans_matrices...']);
			s=strsplit(img_path,filesep);
			dist_trans_matrices={};
			c=1;
			for c=1:length(factor_flag)
				if factor_flag(c)
					data=load(['./temp/dist_trans_matrix/' s{end} '__dist_trans_matrix_' int2str(c) '.mat']);
					dist_trans_matrices{end+1}=data.m;
				end
			end
			%read back centers
			centers=dlmread(centers_filename);
			%make factor matrix
			dist_trans_matrices_num=size(dist_trans_matrices,2);
			centers_num=size(centers,1);
			factor_matrix=zeros(centers_num,dist_trans_matrices_num);
		    for i=1:centers_num
		        for j=1:dist_trans_matrices_num
		            factor_matrix(i,j)=dist_trans_matrices{j}(centers(i,1),centers(i,2),centers(i,3));
		        end
		    end
		    %save factor_matrix into txt file
		    pre_name=strsplit(img_path,filesep);
			factor_matrix_filename = ['./temp/others/' pre_name{end} '__factor_matrix.txt'];
			disp(['saving factor matrix...']);
			save(factor_matrix_filename,'factor_matrix','-ASCII','-double');
		end
	end
end

function save3Dmatrix(m,img_path,dm_n)
	s=strsplit(img_path,filesep);
	filename=['./temp/dist_trans_matrix/' s{end} '__dist_trans_matrix_' int2str(dm_n) '.mat'];
	disp(['saving ' filename]);
	save(filename,'m');
end