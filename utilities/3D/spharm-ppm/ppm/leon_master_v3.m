function leon_master_v3()

imgs_path='/Users/piepi/CellOrganizer/Josh/human*brain*imaging*2nd*round/*-40x*';
% imgs_path='/home/xlu2/Josh/human_brain_multiplex_exm_imaging/human*brain*imaging*1st*round/*40x001.nd2';
flag_vector={{'puncta','factor'},{},{'puncta','model'},{'objects','factor'}};
aspect=[0.161427354511474 0.161427354511474 0.4];
min_obj_size=4;

preprocessing_temp(imgs_path,flag_vector,aspect,min_obj_size);

end