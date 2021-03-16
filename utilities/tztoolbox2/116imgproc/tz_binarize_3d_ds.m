function tz_binarize_3d_ds(rtdir,dbdir,ext,loadthre,imgtype,med)
%TZ_BINARIZE_3D_DS Binarize all the drugimages
%   TZ_BINARIZE_3D_DS(RTDIR,DBDIR,EXT,LOADTHRE,IMGTYPE,MED)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_binarize_all(rtdir,dbdir,ext,loadthre,imgtype,med)
%OVERVIEW:
%  This function binarize all the drugimages and save them as imgtype
%PARAMETERS:
%   rtdir - root directory
%   dbdir - binary images directory
%   ext - resource image extension
%   loadthre - loading threshold
%   imgtype - saving type
%   med - median filtering or not
%DESCRIPTION:
%   Binarize images under DS tree. 
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   12-MAR-2003 Modified TINGZ
%   24-MAR-2003 Modified TINGZ
%   07-JUL-2003 Modified TINGZ
%   09-JUL-2003 Modified TINGZ
%   16-JUL-2003 Modified TINGZ
%   13-DEC-2003 Modified TINGZ

load diroot.mat

prot_list = ml_dir([rtdir '/' DIDIR])
prot_list = prot_list(3:end)

tz_makedir([rtdir '/' dbdir]);

for prot_no = 1:length(prot_list)
    prot_name = prot_list{prot_no};
    prot_fullpath = [rtdir '/' DIDIR '/' prot_name];
    mask_prot_fullpath = [rtdir '/' DMDIR '/' prot_name];
    bin_prot_fullpath = [rtdir '/' dbdir '/' prot_name];
    if(~exist(bin_prot_fullpath,'dir'))
        tz_makedir( bin_prot_fullpath);
    end
    
    drug_list = ml_dir(prot_fullpath);
    drug_list = drug_list(3:end);
    ncondit = length(drug_list);
    %     ndrugs = ncondit/2;
    
    for drug_no = 1:ncondit
        drug_name = drug_list{drug_no};
        drug_fullpath = [prot_fullpath '/' drug_name];
        mask_drug_fullpath = [mask_prot_fullpath '/' drug_name];
        bin_drug_fullpath = [bin_prot_fullpath '/' drug_name];
        if(~exist(bin_drug_fullpath,'dir'))
            tz_makedir( bin_drug_fullpath);
        end
        
        sample_list=ml_dir(drug_fullpath);
        sample_list=sample_list(3:end);
        
        for sample_no = 1:length(sample_list)
            sample_name = sample_list{sample_no};
            nametest=[sample_name '123456'];
            
            if any(nametest(1:7) ~= 'timeser')
                sample_fullpath = [drug_fullpath '/' sample_name];
                mask_fullpath = [mask_drug_fullpath '/' sample_name '.mat'];
                bin_fullpath=[bin_drug_fullpath '/' sample_name]
                if(~exist(bin_fullpath,'dir'))
                    tz_makedir(bin_fullpath);
                    tz_binarize_3d(sample_fullpath,bin_fullpath,mask_fullpath,ext,loadthre,imgtype,med);
                end
                
            end
        end
        
    end
end
