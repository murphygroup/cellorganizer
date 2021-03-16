function feat_vector_ori = real_im_feat_extract(image,imgcent_coordinate,w,cellnum)


imCategory = 'general';
%w = [1;1;1;1;1;1;1];

feat_vector_ori = getfeatvector(image,imgcent_coordinate,imCategory,w,cellnum);

% End of function
