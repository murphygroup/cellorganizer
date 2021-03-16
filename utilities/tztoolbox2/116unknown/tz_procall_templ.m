% A template for processing all images

addpath /home/tingz/matlab/tz_function
addpath /home/tingz/matlab/tz_constant

load shareddir.mat
eval(['addpath ' SHAREDDIR])

load diroot.mat

dipath = [DIROOT DIDIR];
prot_list = ml_dir(dipath);
prot_list = prot_list(3:end);

for prot_no = 1:length(prot_list)
    prot_name = prot_list{prot_no};
    prot_path = [dipath prot_name '/'];

    drug_list = ml_dir(prot_path);
    drug_list = drug_list(3:end);
    ncondit = length(drug_list);
    ndrugs = ncondit/2;

    for drug_no = 1:length(drug_list)
        drug_name = drug_list{drug_no};
        drug_path = [prot_path drug_name '/']

        % Add processing code here.

    end
end


