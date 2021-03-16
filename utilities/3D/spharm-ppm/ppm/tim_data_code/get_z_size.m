function z_size=get_z_size(img_path)
	reader = bfGetReader(img_path);
	omeMeta = reader.getMetadataStore();
	z_size = omeMeta.getPixelsSizeZ(0).getValue();
end