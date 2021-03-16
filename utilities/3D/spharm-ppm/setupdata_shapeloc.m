% return image paths
function [cell_paths,dna_paths,prot_paths,mask_paths] = setupdata_shapeloc(image_root)
cell_lst = ml_ls(image_root);
cell_paths = {};
dna_paths = {};
prot_paths = {};
mask_paths = {};
for i = 1:length(cell_lst)
    cell_paths{i} = [image_root filesep cell_lst{i} filesep 'cell.tif'];
    dna_paths{i} = [image_root filesep cell_lst{i} filesep 'dna.tif'];
    prot_paths{i} = [image_root filesep cell_lst{i} filesep 'prot.tif'];
    mask_paths{i} = [image_root filesep cell_lst{i} filesep 'mask.tif'];
end

end