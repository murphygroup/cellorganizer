function feat_vector = feat_extract2(image,cellmask, dnamask, featType,imgcent_coordinate,numLevels1,numLevels2)
%modified version of feat_extract for microtubule model, that passes in the
%segmented cell image rather than the number, and loading the image from
%disk
%
%segdna is only necessary if featType = 'distHist'

% featType: 'haralick','sliceIntensity','tas','totvolume','histogram','cubIntensity'
% image2 is the part that you want to ignore. this is exclusively for the feature compute - harIgnoreCent


if ~exist('thr','var')
   thr = 6;
end

if ~exist('numLevels1','var')
   numLevels1 = 10;
end
if ~exist('numLevels2','var')
   numLevels2 = 6;
end

if strcmp(featType,'haralick')
% 		feats = ml_texture_withzeros(uint8(image));
        feat_vector = haralickTexture(image, 1);

elseif strcmp(featType,'harkd2')
		feat_vector = haralickTexture(image, 2);

elseif strcmp(featType,'harkd4')
		feat_vector = haralickTexture(image, 4);

elseif strcmp(featType,'radIntensity')
	feat_vector = radInt(image,imgcent_coordinate);

elseif strcmp(featType,'histpropwithcent')
	
	feat_vector = MI_histprops(image(image>0));

elseif strcmp(featType,'histbinwithcent')

	feat_vector = MI_histin(image(image>0));

elseif strcmp(featType,'totint')
	feat_vector = sum(image(:));


elseif strcmp(featType,'edge')
% this works only for 2D . see feat_extract.m in work30/src/ for 3d version
% 	image = double(uint8(image));
    
        
        feat_vector = edgecalc(image, cellmask);
        image = imresize(image,0.5);
        cellmask = imresize(cellmask, 0.5);
        
        feat_vector = [feat_vector,edgecalc(image, cellmask)];
        
        image = imresize(image,0.5);
        cellmask = imresize(cellmask, 0.5);
        
        feat_vector = [feat_vector,edgecalc(image, cellmask)];
        
elseif strcmp(featType,'haralick2')
%        addpath(genpath('/home/jieyuel/lib_backup/HPA_lib/HPA_lib/SLICfeatures/matlab'))
       feat_vector = haralickTexture(sum(image,3), [1,2,4]);


elseif strcmp(featType,'histpropwithcent2')
%        addpath('/home/jieyuel/lib_backup/Microtubules/nocodazole/tarball')
	image2 = image;
	feat_vector = MI_histprops(image2(image2>0));

	image3 = image;
    DSF = 2;
	image3 = imresize(image3, [size(image3,1) size(image3,2)]/DSF);
    feat_vector = [feat_vector, MI_histprops(image3(image3>0))];

	image4 = image;
    DSF = 4;
    image4 = imresize(image4, [size(image4,1) size(image4,2)]/DSF);
    feat_vector = [feat_vector, MI_histprops(image4(image4>0))];


elseif strcmp(featType,'histbinwithcent2')
%        addpath('/home/jieyuel/lib_backup/Microtubules/nocodazole/tarball')
	feat_vector = MI_histbin(image);


elseif strcmp(featType,'histDistLevel2')
	feat_vector = hist_DistLevel2(image,cellmask);


elseif strcmp(featType,'statDistLevel')
	feat_vector = stat_DistLevel(image,cellmask);


elseif strcmp(featType,'distHist')  %%gus
       [l_img] = getDistLevel2_gus2(image,cellmask, dnamask,thr,numLevels2);

       feat_vector = zeros(1,thr);
       for i=1:thr
           feat_vector(i) = nansum(image(l_img == i));
           %feat_vector(i) = nanmean(image(ww));
       end


end % If ends
end % Function ends

