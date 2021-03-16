function leon_master_v2()

imgs_path='/Users/piepi/CellOrganizer/Josh/human*brain*imaging*1st*round/*.nd2';
flag_vector=[18,9,17,18];
aspect=[0.161427354511474 0.161427354511474 0.4];

%check and silence the image that hasn't 4 channels or only one z
imgs_path2=ml_ls(imgs_path);
for i=1:length(imgs_path2)
	img_path=imgs_path2{i};
	reader = bfGetReader(img_path);
	omeMeta = reader.getMetadataStore();
	c_size = omeMeta.getPixelsSizeC(0).getValue();
	z_size = omeMeta.getPixelsSizeZ(0).getValue();
	if (c_size~=4 || z_size<=1) && img_path(end)~='_'
		movefile(img_path,[img_path '_']);
	end
end

factor_matrices=preprocessing(imgs_path,flag_vector,aspect);

end