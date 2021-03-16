function shenj_master(tmp_dir,tiff_dir)
[root_path,diary_name,options,round_num,intermediate_output_path,username,final_output_path,python_path]=ppm_setup(tmp_dir);
% get image paths
% tiff_dir = '/Users/auas/tmp_dir/ometiff';
a = dir([tiff_dir,'/*.tif']);
fn = {a.name};
fd = {a.folder};
imgs_paths = {};
for i=1:length(fn)
    t = [fd{i},'/',fn{i}];
    imgs_paths{i}=t;
end

masks_paths=get_ometiff_mask_paths(imgs_paths,intermediate_output_path);
flag_vector={{'boundaries',{'factor'}},{'nucleus',{'factor'}},{'puncta',{'model'}}}; %./cell1.ome.tif
tic;get_all_distance_transformation_matrix(imgs_paths,masks_paths,intermediate_output_path,flag_vector,options);disp('get_all_distance_transformation_matrix end !');toc;
get_all_models(imgs_paths,masks_paths,intermediate_output_path,final_output_path,python_path,flag_vector,options);disp('get_all_models end !');
% get_all_models(imgs_paths,masks_paths,intermediate_output_path,final_output_path,python_path,flag_vector,options)
end


function [root_path,diary_name,options,round_num,intermediate_output_path,username,final_output_path,python_path]=ppm_setup(tmp_dir)
	[~,username] = system('whoami');
	if strcmpi(strtrim(username), 'xlu2')
		root_path='/home/xlu2/';
	elseif strcmpi(strtrim(username), 'piepi')
		root_path='/Users/piepi/';
    end
    root_path='/Users/auas/GitHub/';
	[current_path, filename, extension] = fileparts( which(mfilename) );
	cd ([root_path 'cellorganizer3/']);
	setup;
	cd(current_path);

	if ~exist( './diaries' )
		system('mkdir ./diaries');
	end
	options.datetime_str=strrep(datestr(datetime),' ','_');
	diary_name=[pwd filesep 'diaries' filesep 'diary_' options.datetime_str '.txt'];
	% disp(['save diary to: ' diary_name]);
	% aspect=[0.161427354511474 0.161427354511474 0.4]; %leon's data
	options.aspect=[0.049,0.049,0.200]; %shenj's data
	options.min_obj_size=4;
    options.mask_inverted_color_flag=0;
	round_num={'probe_a' 'probe_b' 'probe_c' 'probe_d' 'probe_e' 'probe_f' 'probe_g' 'probe_h' 'probe_i' 'probe_j' 'probe_k' 'probe_l' 'probe_m' 'rxn_1'};
% 	intermediate_output_path='./test_ppt';
    intermediate_output_path = [tmp_dir,'/test_ppt'];
	final_output_path='final_result';
	%python_path=[root_path 'anaconda/bin/python'];
    python_path='/Users/auas/anaconda2/bin/python';
	options.dummy_num='50';
	options.rand_num='70000';
	options.cv_mode='rd_roi';
	options.fold='3';
	options.cv_round='1';
	options.debug=1;
end