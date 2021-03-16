function feat_vector = feat_extract(image,image2,featType,imgcent_coordinate,imnum,numLevels1,numLevels2)
% featType: 'haralick','sliceIntensity','tas','totvolume','histogram','cubIntensity'
% image2 is the part that you want to ignore. this is exclusively for the feature compute - harIgnoreCent

if ~exist('numLevels1','var')
   numLevels1 = 10;
end
if ~exist('numLevels2','var')
   numLevels2 = 6;
end


if strcmp(featType,'haralick')
		feat_vector = MI_feat_extract_truncate(image,[1:13]);


elseif strcmp(featType,'radIntensity')
	feat_vector = radInt(image,imgcent_coordinate);


elseif strcmp(featType,'harkd2')
		image = ml_3dXYresize(image,0.5);
        	feat_vector = MI_feat_extract_truncate(image,[1:13]);


elseif strcmp(featType,'harkd4')
		image = ml_3dXYresize(image,0.5);
	        feat_vector = MI_feat_extract_truncate(image,[1:13]);

elseif strcmp(featType,'histpropwithcent')
	image = uint8(image);
	image(find(image==0)) = [];
	feat_vector = MI_histprops(image);

elseif strcmp(featType,'histbinwithcent')
	image = uint8(image);
        image(find(image==0)) = [];
	feat_vector = MI_histbin(image);

elseif strcmp(featType,'totint')
	feat_vector = sum(image(:));


elseif strcmp(featType,'edge')
% this works only for 2D 
	mfs = floor(size(image,3)/2);
	image = double(uint8(image(:,:,mfs)));
        feat_vector = edgecalc(image);
        image = imresize(image,0.5);
        feat_vector = [feat_vector,edgecalc(image)];
        image = imresize(image,0.5);
	feat_vector = [feat_vector,edgecalc(image)];
        
elseif strcmp(featType,'haralick2')
		feat_vector = haralicktexture3d_mvres(image,[1:13]);
		image = ml_3dXYresize(image,0.5);
		feat_vector = [feat_vector, haralicktexture3d_mvres(image,[1:13])];
		image = ml_3dXYresize(image,0.5);
        	feat_vector = [feat_vector, haralicktexture3d_mvres(image,[1:13])];

elseif strcmp(featType,'histpropwithcent2')
	image = uint8(image);
	image2 = image;
	image2(find(image2==0)) = [];
	feat_vector = MI_histprops(image2);
	image3 = image;
	image3 = ml_3dXYresize(image3,0.5);
        image3(find(image3==0)) = [];
        feat_vector = [feat_vector,MI_histprops(image3)];
	image4 = image;
        image4 = ml_3dXYresize(image4,0.25);
        image4(find(image4==0)) = [];
        feat_vector = [feat_vector,MI_histprops(image4)];

elseif strcmp(featType,'histbinwithcent2')
	feat_vector = MI_histbin(image);

elseif strcmp(featType,'histDistLevel2')
       if ismember(imnum,[1,6,14,19,31,52])  %%This is not for curvelet but for the compatibility with curvelet features.
          image = cat(3,image,zeros(size(image,1),size(image,2),16-size(image,3)));
       end
	feat_vector = hist_DistLevel2(image,imnum);
%	image = ml_3dXYresize(image,0.5);
%	feat_vector = [feat_vector, hist_DistLevel2(image,imnum)];
%	image = ml_3dXYresize(image,0.5);
%     	feat_vector = [feat_vector, hist_DistLevel2(image,imnum)];

elseif strcmp(featType,'statDistLevel')
       if ismember(imnum,[1,6,14,19,31,52])  %%This is not for curvelet but for the compatibility with curvelet features.
          image = cat(3,image,zeros(size(image,1),size(image,2),16-size(image,3)));
       end
	feat_vector = stat_DistLevel(image,imnum);

elseif strcmp(featType,'distHist')  %%gus
       if ismember(imnum,[1,6,14,19,31,52])  %%
          image = cat(3,image,zeros(size(image,1),size(image,2),16-size(image,3)));
       end
       [l_img] = getDistLevel2_gus(image,imnum,numLevels2);

       m = max(l_img(:));
       feat_vector = zeros(1,m);
       for i=1:m
           ww = find(l_img == i);
           feat_vector(i) = nansum(image(ww));
           %feat_vector(i) = nanmean(image(ww));
       end

end % If ends
end % Function ends
