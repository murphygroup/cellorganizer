function tz_mat2img_ds(matdir,imgdir,imgtype)
%TZ_MAT2IMG_DS Converts mat files into images (for drug files)
%   TZ_MAT2IMG_DS(MATDIR,IMGDIR,IMGTYPE) converts all mat files in
%   MATDIR into images with format IMGTYPE saves them into IMGDIR.
%   The structure of directory MATDIR is the same as saving drug
%   images.
%
%   See also TZ_MAT2IMG_DIR

%   ??-???-??? Initial write T. Zhao
%   07-JUL-2003 Modified T. Zhao
%   22-MAY-2004 Modified T. Zhao
%       - change old function name
%   Copyright (c) Murphy Lab, Carnegie Mellon University


load diroot.mat

matpath=[DIROOT '/' matdir];
imgpath=[DIROOT '/' imgdir];
tz_make_dir(imgpath);

prot_list=ml_dir(matpath);
prot_list=prot_list(3:end);

for prot_no=1:length(prot_list)
    prot_name=prot_list{prot_no};
    prot_path=[matpath '/' prot_name];
    
    prot_path2=[imgpath '/' prot_name];
    tz_make_dir(prot_path2);
    
    drug_list=ml_dir(prot_path);
    drug_list=drug_list(3:end);
    ncondit=length(drug_list);
    ndrugs=ncondit/2;
    
    for drug_no=1:length(drug_list)
        drug_name=drug_list{drug_no};
        drug_path=[prot_path '/' drug_name]
    
        drug_path2=[prot_path2 '/' drug_name]
        tz_make_dir(drug_path2);
        
        % Add processing code here.
        tz_mat2img_dir(drug_path,drug_path2,imgtype);
        
    end
end

