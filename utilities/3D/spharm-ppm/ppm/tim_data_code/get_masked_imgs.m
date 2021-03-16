function [masked_imgs,ROIIDs]=get_masked_imgs(img,labeled_mask,maskIDs)
	ROIIDs={};
	masked_imgs={};
	for val=1:numel(maskIDs)
	    ROI_mask=labeled_mask;
	    if length(find(ROI_mask==maskIDs(val)))>=100
			ROI_mask(ROI_mask~=maskIDs(val))=0;
			ROI_mask(ROI_mask==maskIDs(val))=1;
			masked_imgs{end+1}=get_maskedImg(img,ROI_mask);
            ROIIDs{end+1}=val;
		end
	end
end

function masked_img=get_maskedImg(img,mask)
	masked_img=zeros(size(img));
	for z=1:size(img,3)
		masked_img(:,:,z)=img(:,:,z).*double(mask);
	end
end