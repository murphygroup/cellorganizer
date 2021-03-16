[root_path,diary_name,options,round_num,intermediate_output_path,username,final_output_path,python_path]=ppm_setup();
% diary(diary_name);
% i=9; %only test for probe_i
i=1;
% for i=1:length(round_num)
	disp(['running in ' round_num{i} '...']);
	[imgs_paths,masks_paths]=get_paths(round_num,i,username,intermediate_output_path);

% 	% flag_vector={{'objects',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'factor'}}};
% 	% flag_vector={{'objects',{'factor'}},{},{},{'puncta',{'model','factor'}},{}}; %1
% 	% flag_vector={{},{'puncta',{'model','factor'}},{},{},{'puncta',{'factor'}}}; %2
% 	% flag_vector={{'objects',{'factor'}},{'puncta',{'model','factor'}},{},{},{'puncta',{'factor'}}}; %3
%	% flag_vector={{'objects',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}}}; %4
% 	% flag_vector=get_flag_vector(imgs_paths);
	
	flag_vector={{'boundaries',{'factor'}},{'nucleus',{'factor'}},{'puncta',{'model'}}}; %./cell1.ome.tif
	tic;get_all_distance_transformation_matrix(imgs_paths,masks_paths,intermediate_output_path,flag_vector,options);disp('get_all_distance_transformation_matrix end !');toc;
	get_all_models(imgs_paths,masks_paths,intermediate_output_path,final_output_path,python_path,flag_vector,options);disp('get_all_models end !');

% % % end
% % print_total_run_time(diary_name);
% % diary off;

function [root_path,diary_name,options,round_num,intermediate_output_path,username,final_output_path,python_path]=ppm_setup()
	[~,username] = system('whoami');
	if strcmpi(strtrim(username), 'xlu2')
		root_path='/home/xlu2/';
	elseif strcmpi(strtrim(username), 'piepi')
		root_path='/Users/piepi/';
	end
	[current_path, filename, extension] = fileparts( which(mfilename) );
	cd ([root_path 'CellOrganizer/cellorganizer3/']);
	setup;
	cd(current_path);

	if ~exist( './diaries' )
		system('mkdir ./diaries');
	end
	options.datetime_str=strrep(datestr(datetime),' ','_');
	diary_name=[pwd filesep 'diaries' filesep 'diary_' options.datetime_str '.txt'];
	% disp(['save diary to: ' diary_name]);
	% aspect=[0.161427354511474 0.161427354511474 0.4]; %leon's data
	options.aspect=[0.0481,0.0481,0.167847]; %tim's data
	options.min_obj_size=4;
	options.mask_inverted_color_flag=0;
	round_num={'probe_a' 'probe_b' 'probe_c' 'probe_d' 'probe_e' 'probe_f' 'probe_g' 'probe_h' 'probe_i' 'probe_j' 'probe_k' 'probe_l' 'probe_m' 'rxn_1'};
	intermediate_output_path='temp_intermediate';
	final_output_path='final_result';
	python_path=[root_path 'anaconda/bin/python'];
	options.dummy_num='50';
	options.rand_num='70000';
	options.cv_mode='rd_roi';
	options.fold='3';
	options.cv_round='1';
	options.debug=1;
end

function [imgs_paths,masks_paths]=get_paths(round_num,i,username,intermediate_output_path)
	if strcmpi(strtrim(username), 'xlu2')
		imgs_paths_=['/home/tmajaira/projects/seema/images/' round_num{i} '/*.ics'];
		mask_paths_=['/home/tmajaira/projects/seema/images/masks/' erase(round_num{i},'probe_') '*mask.tif'];

		imgs_paths={'/home/tmajaira/projects/seema/images/probe_a/6HPI_Probe_A_stk1_decon.ics'};
		masks_paths={'/home/tmajaira/projects/seema/images/masks/a1mask.tif'};
	elseif strcmpi(strtrim(username), 'piepi')
		imgs_paths_=['/Users/piepi/CellOrganizer/Josh/RRA*2/home/tmajaira/projects/seema/images/' round_num{i} '/*.ics'];
		mask_paths_=['/Users/piepi/CellOrganizer/Josh/RRA*2/home/tmajaira/projects/seema/images/masks/' erase(round_num{i},'probe_') '*mask.tif'];

		imgs_paths={'/Users/piepi/CellOrganizer/Josh/RRA 2/home/tmajaira/projects/seema/images/probe_a/6HPI_Probe_A_stk1_decon.ics'};
		masks_paths={'/Users/piepi/CellOrganizer/Josh/RRA 2/home/tmajaira/projects/seema/images/masks/a1mask.tif'};
	end
	% imgs_paths=ml_ls(imgs_paths_);
	
	% masks_paths=ml_ls(mask_paths_);
	
	%%if no mask
	% mask_paths={};

	% if ometiff
	imgs_paths={'./cell1.ome.tif','./cell2.ome.tif','./cell3.ome.tif'};
	masks_paths=get_ometiff_mask_paths(imgs_paths,intermediate_output_path);
end

function flag_vector_preprocessing=get_flag_vector(imgs_paths)
	reader = bfGetReader(imgs_paths{1});
	omeMeta = reader.getMetadataStore();
	c_size = omeMeta.getPixelsSizeC(0).getValue();
	flag_vector_preprocessing={{'objects','factor'}};
	for j=1:(c_size-1)
		flag_vector_preprocessing{end+1}={'puncta','factor'};
	end
end

function print_total_run_time(diary_name)
	filetext = fileread(diary_name);
	C = strsplit(filetext,'\n');
	B=regexp(C,'\d+.*seconds\.$','once','match');
	B=B(find(~cellfun(@isempty,B)));
	times=0;
	for i=1:length(B)
		b=strsplit(B{i},' ');
		times=times+str2num(b{1});
	end
	disp(['total run time:' num2str(times)]);
end

function c_original=trans2original_position(c,flag_vector_preprocessing)
	count=0;
	for i=1:length(flag_vector_preprocessing)
		if ~isempty(flag_vector_preprocessing)
			count=count+1;
			if count==c
				c_original=i;
				break
			end
		end
	end
end