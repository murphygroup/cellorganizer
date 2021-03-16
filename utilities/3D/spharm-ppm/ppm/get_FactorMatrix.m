function get_FactorMatrix(imgs_paths,flag_vector)
%only get factor matrix

	%loop each image
	factor_matrices=[];
	if ~exist( './temp/others' )
		mkdir ./temp/others
	end
	
	for f=1:length(imgs_paths)
	% f=1;
		img_path=imgs_paths{f};
		factor_flag=[];
		%loop each channel
		for c=1:length(flag_vector)
			flag=flag_vector{c};
			if ~isempty(flag)
				if strcmp(flag{2},'model')
					% celldisp(flag);
					factor_flag=[factor_flag,0];
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
                        if strcmp(flag{1},'objects_local')
							[~,centers]=find_objects_local(img,'objects_center',min_obj_size);
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
				if strcmp(flag{2},'factor')
					factor_flag=[factor_flag,1];
				end
			else
				factor_flag=[factor_flag,0];
			end
		end

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