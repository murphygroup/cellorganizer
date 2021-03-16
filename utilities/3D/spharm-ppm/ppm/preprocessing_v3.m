function preprocessing(imgs_paths,flag_vector,aspect,min_obj_size)
%only calculate transformation matrix
	
	%loop each image
	if ~exist( './temp/dist_trans_matrix' )
		system('mkdir -p ./temp/dist_trans_matrix');
	end
	for f=1:length(imgs_paths)
	% f=1;
		img_path=imgs_paths{f};
		dist_trans_matrices={};
		%loop each channel
		for c=1:length(flag_vector)
			flag=flag_vector{c};
			if c==1
				reader = bfGetReader(img_path);
				omeMeta = reader.getMetadataStore();
				z_size = omeMeta.getPixelsSizeZ(0).getValue();
			end
			img=imreadBF(img_path,1:z_size,1,c);
			if strcmp(flag{1},'puncta')
				[puncta_img,~]=find_punta(img,'puncta_img');
				disp([img_path '--channel:' c ' transformation matrix caculation start...']);
				dist_trans_matrices{end+1}=bwdistsc(puncta_img,aspect);
				save3Dmatrix(dist_trans_matrices{end},img_path,c);
			end
			if strcmp(flag{1},'objects')
				[objects_img,~]=find_objects(img,'objects_img',min_obj_size);
				disp([img_path '--channel:' c ' transformation matrix caculation start...']);
				dist_trans_matrices{end+1}=bwdistsc(objects_img,aspect);
				save3Dmatrix(dist_trans_matrices{end},img_path,c);
			end
		end
	end
end

function save3Dmatrix(m,img_path,dm_n)
	s=strsplit(img_path,filesep);
	filename=['./temp/dist_trans_matrix/' s{end} '__dist_trans_matrix_' int2str(dm_n) '.mat'];
	disp(['saving ' filename]);
	save(filename,'m');
end