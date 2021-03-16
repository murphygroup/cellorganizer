function master()

%read image
current_path = which(mfilename);
[current_path, filename, extension] = fileparts( current_path );
cd(current_path);

% imgs_path='/Users/piepi/CellOrganizer/Josh/human brain imaging 1st round/AR3-Trs pH 10-40x001.nd2';
imgs_path='/home/xlu2/Josh/human_brain_multiplex_exm_imaging/human brain imaging 1st round/AR3-Trs pH 10-40x001.nd2';
reader = bfGetReader(imgs_path);
omeMeta = reader.getMetadataStore();
% c_size = omeMeta.getPixelsSizeC(0).getValue();
z_size = omeMeta.getPixelsSizeZ(0).getValue();
model_channel_imgs=imreadBF(imgs_path,1:z_size,1,2); %[vol]=imreadBF(datname,zplanes,tframes,channel)
other_channel_imgs={};
other_channel_imgs{1}=imreadBF(imgs_path,1:z_size,1,1);
other_channel_imgs{2}=imreadBF(imgs_path,1:z_size,1,3);
other_channel_imgs{3}=imreadBF(imgs_path,1:z_size,1,4);
puncta_or_threshold_flags=[0,1,0];

factor_matrix=calculate_factor(model_channel_imgs,other_channel_imgs,puncta_or_threshold_flags);

end

%imshow(dist_trans_matrices3.temp_dist_trans_matrices(:,:,16)) %%%check object