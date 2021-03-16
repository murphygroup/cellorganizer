function feats=tz_objdnadist(objects,dnaobj)
%TZ_OBJDNADIST Unknown.
%   FEATS = TZ_OBJDNADIST(OBJECTS,DNAOBJ)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function feats=tz_objfeat_3d(object)
%
%OVERVIEW:
%   calculate features for 3d objects
%

dnaobj=double(dnaobj);
dnacof=sum([dnaobj(:,1).*dnaobj(:,4),dnaobj(:,2).*dnaobj(:,4),dnaobj(:,3).*dnaobj(:,4)],1)/sum(dnaobj(:,4));

for i=1:length(objects)
    
    object=double([objects{i}.voxels;objects{i}.gray]');
    
    objcof(i,:)=sum([object(:,1).*object(:,4),object(:,2).*object(:,4),object(:,3).*object(:,4)],1)/sum(object(:,4));
    
    %3D-SOF2.12: Horizontal 2D Euclidean distance between object COF and DNA COF
    feats(i)=sqrt(sum((objcof(i,1:2)-dnacof(1:2)).^2));  
   
end
feats=feats';

