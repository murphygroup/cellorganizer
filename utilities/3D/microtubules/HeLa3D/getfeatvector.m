function [feat_vector_ori,idxes] = getfeatvector(im4,imgcent_coordinate,imCategory,w,cellnum)
% Get any images feature vector

feat_vector_ori = [];
im5 = [];

if (w(1) == 1)
        idxes{1}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'haralick',imgcent_coordinate)];
        idxes{1}(2) = size(feat_vector_ori,2);
end


if (w(2) == 1)
        idxes{2}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'histpropwithcent',imgcent_coordinate)];
        idxes{2}(2) = size(feat_vector_ori,2);
end

if (w(3) == 1)
        idxes{3}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'harkd2',imgcent_coordinate)];
        idxes{3}(2) = size(feat_vector_ori,2);
end

if (w(4) == 1)
        idxes{4}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'harkd4',imgcent_coordinate)];
        idxes{4}(2) = size(feat_vector_ori,2);
end

if (w(5) == 1)
        idxes{5}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'totint',imgcent_coordinate)];
         idxes{5}(2) = size(feat_vector_ori,2);
end

if (w(6) == 1)
        idxes{6}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori = [feat_vector_ori,feat_extract(im4,[],'radIntensity',imgcent_coordinate)];
        idxes{6}(2) = size(feat_vector_ori,2);
end

if (w(7) == 1)
        idxes{7}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,[],'edge',imgcent_coordinate)];
        idxes{7}(2) = size(feat_vector_ori,2);
end

if (w(8) == 1)
	idxes{8}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'haralick2',imgcent_coordinate)];
	idxes{8}(2) = size(feat_vector_ori,2);
end

if (w(9) == 1)
	idxes{9}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'histpropwithcent2',imgcent_coordinate)];
	idxes{9}(2) = size(feat_vector_ori,2);
end

if (w(10) == 1)
	idxes{10}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'histbinwithcent2',imgcent_coordinate)];
	idxes{10}(2) = size(feat_vector_ori,2);
end

if (w(11) == 1)
	idxes{11}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'histDistLevel2',imgcent_coordinate,cellnum)];
	idxes{11}(2) = size(feat_vector_ori,2);
end

if (w(12) == 1)
	idxes{12}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'statDistLevel',imgcent_coordinate,cellnum)];
	idxes{12}(2) = size(feat_vector_ori,2);
end

if (w(13) == 1)
	idxes{13}(1) = size(feat_vector_ori,2)+1;
        feat_vector_ori  = [feat_vector_ori,feat_extract(im4,im5,'distHist',imgcent_coordinate,cellnum)];
	idxes{13}(2) = size(feat_vector_ori,2);
end

end % function ends
