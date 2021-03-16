function [big_matrix, name_matrix] = tz_loadfeats_ds(dfpath)
%TZ_LOADFEATS_DS Load features for drug images.
%   BIG_MATRIX = TZ_LOADFEATS_DS(DFPATH) loads SLF from the directory
%   DFPATH. BIGMATRIX is a 3D cell array of feature matrix.
%
%   [BIG_MATRIX,NAME_MATRIX] = TZ_LOADFEATS_DS(DFPATH) also returns
%   group name for each feature matrix.

%   07-JUN-2003 Initial write TINGZ
%   08-JUN-2003 Modified TINGZ
%   08-JUN-2003 Modified TINGZ
%   14-DEC-2003 Modified TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

prot_list=ml_dir(dfpath);
prot_list=prot_list(3:end);

big_matrix={};
name_matrix={};

for prot_no=1:length(prot_list)
    prot_name=prot_list{prot_no};
    prot_path=[dfpath '/' prot_name];
    
    drug_list=ml_dir(prot_path);
    drug_list=drug_list(3:end);
    ncondit=length(drug_list);
    ndrugs=ncondit/2;
    
    for drug_no=1:length(drug_list)
        drug_name=drug_list{drug_no};
        drug_path=[prot_path '/' drug_name]
        
        % Add processing code here.
        sample_list=tz_readdir([drug_path '/*.mat']);
        
        features=[];
        
        for sample_no=1:length(sample_list)
            load(sample_list{sample_no});
            features=[features;image_features];
        end
        
        if (drug_no > ndrugs)
            big_matrix{1, drug_no-ndrugs, prot_no} = features;
            name_matrix{1, drug_no-ndrugs, prot_no} = {prot_name, drug_name};
        else
            big_matrix{2, drug_no, prot_no} = features;
            name_matrix{2, drug_no, prot_no} = {prot_name, drug_name};
        end     
    end
end

