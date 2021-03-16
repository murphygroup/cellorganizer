function [ofeat,oname] = features_overlap( prot, procprot, procref, refname)
% Takes uint8 images, already thresholded but with intensity

if ~exist('refname','var')
    refname = 'ref';
end

% Threshold by Otsu's method

obj_prot = procprot>0;
obj_ref = procref>0;

area_prot = sum(obj_prot(:));
area_ref = sum(obj_ref(:));

overlap = obj_prot & obj_ref;
areaOverlap = sum(overlap(:));

intersect_area_to_prot_area = areaOverlap / area_prot;
intersect_area_to_ref_area = areaOverlap / area_ref;

intersect_inten_to_prot_inten = sum(prot(obj_ref)) / sum(prot(:));

dist_ref = bwdist(obj_ref);
prot_nuc_pixel_distance_to_prot_area = sum(dist_ref(obj_prot))/area_prot;


ofeat = [intersect_area_to_prot_area intersect_area_to_ref_area ...
    intersect_inten_to_prot_inten ...
    prot_nuc_pixel_distance_to_prot_area];
    
oname = {'Fraction of protein object area that overlaps the reference area', ...
         'Fraction of referece area that overlaps the protein area', ...
         'Fraction of protein intensity that overlaps the reference area', ...
         'Avg distace of protein containing pixel to nearest pixel of reference'};

for i=1:length(oname)
    oname{i} = [refname ':' oname{i}];
end


return

