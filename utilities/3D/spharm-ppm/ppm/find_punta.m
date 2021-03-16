function [puncta_img,centers]=find_punta(img,flag)
	disp(['finding puncta...']);

	[fg_denoise] = seg(img,1.5,5);
	objs = bwconncomp(fg_denoise);
	centers = struct2cell(regionprops(objs,'Centroid'));
	centers = cellfun(@int16,centers,'UniformOutput',false);

	if strcmp(flag,'puncta')
		puncta_img = [];
	end

	if strcmp(flag,'puncta_img')
	    puncta_img=zeros(size(img));
	    for i=1:size(centers,2)
	        puncta_img(centers{i}(1),centers{i}(2),centers{i}(3))=1;
	    end
	end
	
    disp(['puncta finding complete']);
end

function [fg_denoise] = seg(improt,f1,f2)
	fg_denoise = zeros(size(improt));
	for z = 1:size(improt,3)
	    curslice = double(improt(:,:,z));

	    f = fspecial('gaussian', 15, f1);
	    improt_f = imfilter(curslice,f);
	    f = fspecial('gaussian', 15, f2);
	    im_f = imfilter(curslice, f);

	    improt_fg = improt_f - im_f;
	    improt_fg(improt_fg < 0) = 0;

	    improt_fg_denoise = improt_fg;
	    improt_fg_denoise(improt_fg_denoise < ml_rcthreshold(uint8(improt_fg_denoise))) = 0;
	    improt_fg_denoise = improt_fg_denoise.* bwmorph(improt_fg_denoise, 'clean', inf);

	    fg_denoise(:,:,z) = improt_fg_denoise;
	end
end