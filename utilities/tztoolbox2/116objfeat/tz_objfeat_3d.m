function feats=tz_objfeat_3d(object,dnaobj,res)
%TZ_OBJFEAT_3D Calculate 3D object level features. 
%

%function feats=tz_objfeat_3d(object)
%
%OVERVIEW:
%   calculate features for 3d objects
%PARAMETERS:
%   object - objects found by binarizingu
%   dnaobj - dna objects
%   imgsize - the original image size
%RETURN:
%   feats - features
%DESCRIPTION:
%   features:
%       3D-SOF1.1: Number of voxels in object
%       3D-SOF1.2: 3D Euclidean distance between obj COF and DNA COF
%       3D-SOF1.3: Fraction of object voxels overlapping with DNA
%       3D-SOF1.4: A measure of eccentricity of the object
%       3D-SOF1.5: Number of holes in the object
%       3D-SOF1.6: A measure of roundness of the object
%       3D-SOF1.7: The length of the object's skeleton by homotopic thinning
%       3D-SOF1.8: The ratio of skeleton length to the area of the convex hull of the skeleton
%       3D-SOF1.9: The fraction of object pixel contained within the skeleton
%       3D-SOF1.10: The fraction of object fluorescence contained within the skeleton
%       3D-SOF1.11: The ratio of the number of branch points in skeleton to length of skeleton
%       3D-SOF2.12: Horizontal 2D Euclidean distance between object COF and DNA COF
%       3D-SOF2.13: Vertical 1D signal distance between object COF and DNA COF
%HISTORY:
%   27-NOV-2004 Initial write TINGZ
%

error(tz_genmsg('of','tz_objfeat_3d'));

% dnaobj=double(dnaobj);
dnacof=sum([dnaobj(:,1).*dnaobj(:,4),dnaobj(:,2).*dnaobj(:,4),dnaobj(:,3).*dnaobj(:,4)],1)/sum(dnaobj(:,4));
dnaind=sub2ind(imgsize,dnaobj(:,1),dnaobj(:,2),dnaobj(:,3));
dnaimg=tz_obj2img(dnaobj,imgsize,{'3d','bn'});
emptyimg=zeros(imgsize);

for i=1:length(objects)
    
%     object=[double(objects{i}.voxels);double(objects{i}.gray)]';
%     objind=sub2ind(imgsize,object(:,1),object(:,2),object(:,3));
%        
%     
%     %3D-SOF1.1: Number of voxels in object
%     feats(i,1)=size(object,1);
%     
%     %3D-SOF1.2: 3D Euclidean distance between obj COF and DNA COF
%     %object cof
%     objcof(i,:)=sum([object(:,1).*object(:,4),object(:,2).*object(:,4),object(:,3).*object(:,4)],1)/sum(object(:,4));
%     feats(i,2)=0; 
%     
%     
%     objimg=emptyimg;
%     objimg(objind)=objects{i}.gray;
%     %bobjimg=objimg>0;
% %     tic
%     %3D-SOF1.3: Fraction of object voxels overlapping with DNA
%     feats(i,3)=sum(dnaimg(objind));
%     
% %     toc
%     %feats(i,3)=tz_countsame(dnaind,objind)/feats(i,1);
%     
%    % objmask=sum(objimg,3);
%     
% %     tic
%     %remove redundant background
%     sobjimg=objimg(:,:,min(object(:,3)):max(object(:,3)));
%     sobjimg=sobjimg(min(object(:,1)):max(object(:,1)),:,:);
%     sobjimg=sobjimg(:,min(object(:,2)):max(object(:,2)),:);
%     %xproj=sum(objmask,1);
%     %yproj=sum(objmask,2);
% %     toc
%     %zproj=sum(sum(objimg,1),2);
% 
% %     sobjmask(yproj==0,:)=[];
% %     sobjmask(:,xproj==0)=[];
%     sobjmask=sum(sobjimg,3);
%     
% %     tic
%     %3D-SOF1.4: A measure of eccentricity of the object
%     obj_feature = imfeature(sobjmask>0,'Eccentricity');
%     feats(i,4)=obj_feature.Eccentricity;
% %     toc
%     
%     %3D-SOF1.5: Number of holes in the object
%     feats(i,5)=objects{i}.n_holes;
% %     tic
%     %3D-SOF1.6: A measure of roundness of the object
%     objimageperim = bwarea(bwperim(sobjmask>0)) ;
%     feats(i,6) = (objimageperim^2)/(4*pi*feats(i,1)) ;
% %     toc
%     
% %     tic
%     %3D-SOF1.7: The length of the object's skeleton by homotopic thinning
%     %3D-SOF1.8: The ratio of skeleton length to the area of the convex hull of the skeleton
%     %3D-SOF1.9: The fraction of object pixel contained within the skeleton
%     %3D-SOF1.10: The fraction of object fluorescence contained within the skeleton
%     %3D-SOF1.11: The ratio of the number of branch points in skeleton to length of skeleton
%     zproj=sum(sum(sobjimg,1),2);
%     [maxinten,sel]=max(zproj);
%     feats(i,7:11)=tz_objskelfeats(sobjimg(:,:,sel(1)));    
% %     toc
%     
%     %3D-SOF2.12: Horizontal 2D Euclidean distance between object COF and DNA COF
%     feats(i,12)=sqrt(sum((objcof(i,1:2)-dnacof(1:2)).^2));  
%     
%     %3D-SOF2.13: Vertical 1D signal distance between object COF and DNA COF
%     feats(i,13)=objcof(i,3)-dnacof(3);
    feats(i,1)=sum(object(:,4));
end
% tic
% feats(:,2)=mv_eucdist(dnacof',objcof')';
% toc

