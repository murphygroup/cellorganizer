function [puncta_img,centers]=find_objects(img,flag,min_obj_size)
	disp(['finding objects...']);

	BW=zeros(size(img));
	for i=1:size(img,3)
		BW(:,:,i)=imbinarize(int16(img(:,:,i)));
	end
	
	objs = bwconncomp(BW);
	if strcmp(flag,'objects_center')
		centers_all = struct2cell(regionprops(objs,'Centroid'));
		centers_all = cellfun(@int16,centers_all,'UniformOutput',false);
		numPixels = cellfun(@numel,objs.PixelIdxList);
        idx=find(numPixels>min_obj_size);
		centers=centers_all(idx);
		puncta_img = [];
	end

	if strcmp(flag,'objects_img')
        numPixels = cellfun(@numel,objs.PixelIdxList);
        idx=find(numPixels<=min_obj_size);
        for i=1:length(idx)
        	BW(objs.PixelIdxList{idx(i)}) = 0;
        end
        puncta_img=BW;
        centers={};
	end

	disp(['objects finding complete']);
end