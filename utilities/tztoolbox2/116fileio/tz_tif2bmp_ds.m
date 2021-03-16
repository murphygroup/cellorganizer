function tz_tif2bmp_ds(rtdir,dbdir,ext,is_3d,imgtype)
%TZ_TIF2BMP_DS Obsolete

%function tz_tif2bmp_ds(rtdir,dbdir,ext,is_3d,imgtype)
%
%OVERVIEW
%   This function convert tif images into images of another format for ds
%PRAMETERS:
%   rtdir - root directory
%   dbdir - destiantion directory
%RETURN
%
%DESCRIPTION
%
%HISTORY
%   ??-???-???? Initial write TINGZ
%   12-JUL-2003 Modified TINGZ
%   14-DEC-2003 Modified TINGZ

load diroot.mat

prot_list = ml_dir([rtdir '/' DIDIR])
prot_list = prot_list(3:end)

tz_make_dir([rtdir '/' dbdir]);

for prot_no = 1:length(prot_list)
    prot_name = prot_list{prot_no};
    prot_fullpath = [rtdir '/' DIDIR '/' prot_name];
    mask_prot_fullpath = [rtdir '/' DMDIR '/' prot_name];
    bmp_prot_fullpath = [rtdir '/' dbdir '/' prot_name];
    if(~exist(bmp_prot_fullpath,'dir'))
        tz_make_dir( bmp_prot_fullpath);
    end
    
    drug_list = ml_dir(prot_fullpath);
    drug_list = drug_list(3:end);
    ncondit = length(drug_list);
    %     ndrugs = ncondit/2;
    
    for drug_no = 1:ncondit
        drug_name = drug_list{drug_no};
        drug_fullpath = [prot_fullpath '/' drug_name];
        mask_drug_fullpath = [mask_prot_fullpath '/' drug_name];
        bmp_drug_fullpath = [bmp_prot_fullpath '/' drug_name];
        if(~exist(bmp_drug_fullpath,'dir'))
            tz_make_dir( bmp_drug_fullpath);
        end
        
        sample_list=ml_dir(drug_fullpath);
        sample_list=sample_list(3:end);
        
        for sample_no = 1:length(sample_list)
            sample_name = sample_list{sample_no};
            nametest=[sample_name '123456'];
            
            if any(nametest(1:7) ~= 'timeser')
                sample_fullpath = [drug_fullpath '/' sample_name]
                mask_fullpath = [mask_drug_fullpath '/' sample_name '.mat'];
                bmp_fullpath=[bmp_drug_fullpath '/' sample_name]
                if(~exist(bmp_fullpath,'dir'))
                    tz_make_dir(bmp_fullpath);
                    tz_tif2bmp(sample_fullpath,bmp_fullpath,ext,is_3d,imgtype);
                end
                
            end
        end
        
    end
end