function [puncta_img,centers]=find_objects(img,flag,min_obj_size)
	disp(['finding objects...']);
    flag_trival=0;
    flag_binarized=0;
    BW=imbinarize(img);
	% BW=zeros(size(img));
	% for i=1:size(img,3)
	% 	BW(:,:,i)=imbinarize(int16(img(:,:,i)));
	% 	if length(unique(BW(:,:,i)))==1 && length(unique(img(:,:,i)))>2
	% 		if flag_trival==0
	% 			warning('Current image cannot be binarized using imbinarize function, try to use a trival method instead.');
	% 		end
	% 		img2=img(:,:,i);
	% 		img2(find(img2~=0))=1;
	% 		BW(:,:,i)=int16(img2);
 %            flag_trival=1;
 %        end
 %        if length(unique(BW(:,:,i)))==1 && length(unique(img(:,:,i)))==2
	% 		if flag_binarized==0
	% 			warning('Current image has been a binarized image, do not need to binarize.');
 %            end
	% 		BW(:,:,i)=int16(img(:,:,i));
 %            flag_binarized=1;
	% 	end
	% end
	
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