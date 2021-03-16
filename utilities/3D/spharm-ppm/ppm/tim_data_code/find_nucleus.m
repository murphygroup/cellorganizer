function [nucleus_img,centers]=find_nucleus(img,flag,min_obj_size)
	disp(['finding nucleus...']);
    flag_trival=0;
    flag_binarized=0;
    BW=imbinarize(img,graythresh(img)*50);
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
	if strcmp(flag,'nucleus_center')
		centers_all = struct2cell(regionprops(objs,'Centroid'));
		centers_all = cellfun(@int16,centers_all,'UniformOutput',false);
		numPixels = cellfun(@numel,objs.PixelIdxList);
        idx=find(numPixels>min_obj_size);
		centers=centers_all(idx);
		nucleus_img = [];
	end

	if strcmp(flag,'nucleus_img')
        numPixels = cellfun(@numel,objs.PixelIdxList);
        idx=find(numPixels<=min_obj_size);
        for i=1:length(idx)
        	BW(objs.PixelIdxList{idx(i)}) = 0;
        end
        nucleus_img=BW;
        centers={};
    end
    % imshow(nucleus_img(:,:,9))
	disp(['nucleus finding complete']);
end