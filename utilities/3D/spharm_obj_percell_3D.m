function cellparam = spharm_obj_percell_3D(imdna_path,imcell_path,improt_path,immask_path,options,seg)
% segment 3D protein image and prepare for
% building spharm model and ppm model
% 12/5/2020 R.F. Murphy - consolidate calls to find_objects_local and add max_obj_size
% 1/26/2021 R.F. Murphy - save obj centers; calculate and save normalized distance to nucleus
% 1/27/2021 R.F. Murphy - also save normdists and angles; don't suppress
%                           objects with zero distcodes
% 2/1/2021 R.F. Murphy - mask protein image before object finding

min_obj_size = options.min_obj_size;
max_obj_size = options.max_obj_size;
sigma = options.local_thresholding_sigma;
thresPerc = options.object_detection_thresPerc;
s = options.protein_image_path;
t = split(s,'/');
%hard code for now
if strcmp(t{length(t)},'prot.tif')
    cellName = t{length(t)-1};
else
    cellName = t{length(t)};
end

% load image - only need the protein image
%imdna = loadImage(imdna_path, [1,1,1]);
%imcell = loadImage(imcell_path, [1,1,1]);
improt = loadImage(improt_path, [1,1,1]);
%immask = loadImage(immask_path, [1,1,1]);
improt(~seg.cell) = 0;
[obj_img,numPixels,centers] = find_objects_local(improt,'objects_img',min_obj_size,max_obj_size,sigma,thresPerc);

ppm_dir = [pwd filesep 'ppm_input']; % replace improt image with puncta as center
spharm_dir = [pwd filesep 'spharm_input']; % save all the objects as 3D tif file
[status, msg, msgID] = mkdir(ppm_dir);
[status, msg, msgID] = mkdir(spharm_dir);
save_spharm_obj(obj_img,min_obj_size,cellName) % save all the objects into 3D tif files
% save ometiff files for objects
% generate puncta image for image:

%Serena - 08/20
%cellparam.seg.keep_centers=keep_centers;
cellparam.seg.obj_img=obj_img;
cellparam.seg.centers=centers;
cellparam.seg.sizes=numPixels;
% 1/27/21 RFM this code not needed because calculated below including angles
%celldist = bwdist(imcell);
%nucdist = bwdist(imdna);
%for i=1:length(centers)
%    cdist = celldist(centers{i});
%    ndist = nucdist(centers{i});
%    normdists{i} = ndist/(cdist+ndist);
%end
%cellparam.seg.normdists=normdists;


%TO FIND BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = ml_initparam(options, struct('verbose', false, 'debug', false));

beta = [];

options = ml_initparam(options, struct('use_geodesic_distance', false));


%Serena 03/21
se = strel('line', 3,5);
seg.cell=imfill(seg.cell,'holes');
cellerode = imerode(seg.cell, se);
cellEdge = logical(abs(imsubtract(cellerode, seg.cell)));
cellEdge=imfill(cellEdge,'holes');


seg.nuc=imfill(seg.nuc,'holes');
nucerode=imerode(seg.nuc,se);
nucEdge= logical(abs(imsubtract(nucerode, seg.nuc)));

%nucEdge = bwperim(seg.nuc,18);
%cellEdge = bwperim(seg.cell,18);
obj_center_image = zeros(size(cellEdge));
            
for p=1:size(centers)
    objcenters=centers(p,:);
    for t = 1:size(objcenters(:))
         obj_center_image(objcenters{1,t}(1), objcenters{1,t}(2), objcenters{1,t}(3)) = 1;
    end
end

disp('Find normalized coordinates for positions of objects');
[distcodes,coords,angles] = ...
      ml_celldistcode(cellEdge,nucEdge,obj_center_image,1,'all', options.use_geodesic_distance);

% don't delete these here because it will throw off the indexing of images
% after shape parameterization
%nullidx = find(distcodes(:,1)==0 & distcodes(:,2)==0);
%distcodes(nullidx,:)=[];
%angles.theta(nullidx,:)=[];
%angles.phi(nullidx,:)=[];
 
%normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);

for i=1:length(distcodes(:,1))
    %Serena 03/21
    if and(distcodes(i,1)==0,distcodes(i,2)==0)
        normdists(i,1)=0;
        continue;
    end
    normdists(i,1)=distcodes(i,1)/sum(abs(distcodes(i,1:2)),2);
end

disp('Mapping positions')
x = ml_mappos([normdists angles.theta angles.phi]);

cellparam.seg.distcodes = distcodes;
cellparam.seg.coords = coords;
cellparam.seg.angles.theta = angles.theta;
cellparam.seg.angles.phi = angles.phi;
cellparam.seg.normdists = normdists;
cellparam.seg.mappos_x = x;
% don't fit for each cell, do the fit later for all cells
%if size(x,1) > size(x,2)
%    beta(size(beta,1)+1,:) = ml_logreg(x,distcodes(:,3))';
%else
%    beta = [];
%end
%cellparam.seg.beta=beta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cellparam.msg = ['done preprocessing spharm objects'];

end

