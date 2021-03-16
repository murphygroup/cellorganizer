function [labeled_mask,maskIDs]=get_labeled_mask(mask_path,mask_inverted_color_flag)
	mask=ml_readimage(mask_path);
	if mask_inverted_color_flag
		mask=imcomplement(mask);
	end
	labeled_mask=bwlabeln(mask);
	maskIDs=unique(labeled_mask);
	maskIDs=maskIDs(maskIDs>0);
end