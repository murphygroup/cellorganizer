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
datetime_str=strrep(datestr(datetime),' ','_');
diary_name=[pwd filesep 'diaries' filesep 'diary_' datetime_str '.txt'];
disp(['save diary to: ' diary_name]);
diary(diary_name);
% aspect=[0.161427354511474 0.161427354511474 0.4]; %leon's data
aspect=[0.0481,0.0481,0.167847]; %tim's data
min_obj_size=4;
round_num={'probe_a' 'probe_b' 'probe_c' 'probe_d' 'probe_e' 'probe_f' 'probe_g' 'probe_h' 'probe_i' 'probe_j' 'probe_k' 'probe_l' 'probe_m' 'rxn_1'};

i=9; %only test for probe_i
% i=2;
% for i=1:length(round_num)
	disp(['running in ' round_num{i} '...']);
	if strcmpi(strtrim(username), 'xlu2')
		imgs_paths_=['/home/tmajaira/projects/seema/images/' round_num{i} '/*.ics'];
		mask_paths_=['/home/tmajaira/projects/seema/images/masks/' erase(round_num{i},'probe_') '*mask.tif'];
	elseif strcmpi(strtrim(username), 'piepi')
		imgs_paths_=['/Users/piepi/CellOrganizer/Josh/RRA*2/home/tmajaira/projects/seema/images/' round_num{i} '/*.ics'];
		mask_paths_=['/Users/piepi/CellOrganizer/Josh/RRA*2/home/tmajaira/projects/seema/images/masks/' erase(round_num{i},'probe_') '*mask.tif'];
	end
	imgs_paths=ml_ls(imgs_paths_);
	mask_paths=ml_ls(mask_paths_);
	%%if no mask
	% mask_paths={};

	%preprocessing
	reader = bfGetReader(imgs_paths{1});
	omeMeta = reader.getMetadataStore();
	c_size = omeMeta.getPixelsSizeC(0).getValue();
	flag_vector_preprocessing={{'objects','factor'}};
	for j=1:(c_size-1)
		flag_vector_preprocessing{end+1}={'puncta','factor'};
	end
	% flag_vector_preprocessing={{},{'puncta','factor'},{'puncta','factor'},{'puncta','factor'},{'puncta','factor'}};
	% tic;preprocessing_tim_v1(imgs_paths,flag_vector_preprocessing,aspect,min_obj_size,mask_paths);disp('preprocessing end !');toc;
	
	%%processing all possible flags
	puncta_positions=get_puncta_positions(flag_vector_preprocessing); %The positions are the potential model channel positions after remove all {} (empty flags)
	for p_ind=1:length(puncta_positions)
		p=puncta_positions(p_ind);
		flag_vectors=get_flag_vectors(flag_vector_preprocessing,p);
		for ii=1:length(flag_vectors)
			% ii=length(flag_vectors); %%% need to comment out
			% ii=length(flag_vectors)-1; %%% need to comment out
			% flag_vectors{ii}={{},{'puncta','factor'},{'puncta','factor'},{},{'puncta','model'}}; %% this is the one cause the logl>0 issue 
			fprintf('\n%i,%i\n',p,ii);
			processing(imgs_paths,flag_vectors{ii},ii,round_num,i,root_path,mask_paths);
			% break
		end
		% break
		print_best_model(flag_vectors,p,datetime_str,i,flag_vector_preprocessing);
		% break
	end
% end
print_total_run_time(diary_name);
diary off;

function processing(imgs_paths,flag_vector,ii,round_num,i,root_path,mask_paths)
%get factor matrix and build model
	celldisp(flag_vector);
	if ii==1
		if exist( './full_model.txt' )
		    delete( './full_model.txt' );
		end
		if exist( './likelihood.txt' )
		    delete( './likelihood.txt' );
		end
	end
	save_factor_channel_num(flag_vector);
	
	%get factor matrix
	tic;get_FactorMatrix_tim(imgs_paths,flag_vector,mask_paths);disp('get factor matrix end !');toc;
	
	%build model
	[~,username] = system('whoami');
	if strcmpi(strtrim(username), 'xlu2')
		folder_name=['/home/tmajaira/projects/seema/images/' round_num{i} '/'];
	elseif strcmpi(strtrim(username), 'piepi')
		% folder_name='/Users/piepi/CellOrganizer/Josh/RRA 2/home/tmajaira/projects/seema/images/probe_i/';
		folder_name=['/Users/piepi/CellOrganizer/Josh/RRA 2/home/tmajaira/projects/seema/images/' round_num{i} '/'];
	end
	regex_input='^[^\.].*.ics$';
	temp_path='CellOrganizer/cellorganizer3/demos/3D/demo3Dpoint_process_module/tim_data_code/temp/';

	tic;
	[status,cmdout] = system([root_path 'anaconda/bin/python model_building_tim_v1.py "' folder_name '" "' regex_input '" "' root_path '" "' temp_path '"'],'-echo');
	if ~isempty(strfind(cmdout,'Traceback (most recent call last):'))
		toc;
		error('Encounter error(s) when build the model !');
	else
		disp('model building end !');
		toc;
	end
end

function print_best_model(flag_vectors,c,datetime_str,i,flag_vector_preprocessing)
	c_original=trans2original_position(c,flag_vector_preprocessing);
	diary_name=[pwd filesep 'diaries' filesep 'result_' datetime_str '_round' num2str(i) '_channel' num2str(c_original) '.txt'];
	disp(['save report to: ' diary_name]);
	diary(diary_name);

	numLines = length(flag_vectors);
	txt_filename='./full_model.txt';
	full_models=readTxt2cell(txt_filename,numLines);

	fid = fopen('likelihood.txt');
	likelihoods = textscan(fid,'%s%s');
	fclose(fid);
	avg=[];
	for j=1:7
		avg=[avg,str2num(likelihoods{2}{2*j-1})];
	end
	best_ind=find(avg==max(avg));

	disp('best model flag vector:');
	celldisp(flag_vectors(best_ind));
	disp(['coefficients: ' full_models{best_ind}]);
	disp(['likelihood: ' likelihoods{2}{2*best_ind-1}]);
	disp(['likelihood sd: ' likelihoods{2}{2*best_ind}]);

	diary off;
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

function flag_vector=insertMC2cellArray(factor_vector,p,flag_vector_preprocessing);
	%before consider empty flag, add MC (model channel flag)
	len=length(factor_vector)+1;
	if p==1
		flag_vector_before{1}={'puncta','model'};
		flag_vector_before(2:len)=factor_vector;
	elseif p==len
		flag_vector_before{p}={'puncta','model'};
		flag_vector_before(1:p-1)=factor_vector;
	else
		flag_vector_before{p}={'puncta','model'};
		flag_vector_before(1:p-1)=factor_vector(1:p-1);
		flag_vector_before(p+1:len)=factor_vector(p:end);
	end
	%add empty flag
	empty_flag=0;
	for i=1:length(flag_vector_preprocessing)
		if isempty(flag_vector_preprocessing{i})
			empty_flag=1;
			len=length(flag_vector_before)+1;
			if i==1
				flag_vector{1}={};
				flag_vector(2:len)=flag_vector_before;
			elseif i==len
				flag_vector{i}={};
				flag_vector(1:i-1)=flag_vector_before;
			else
				flag_vector{i}={};
				flag_vector(1:i-1)=flag_vector_before(1:i-1);
				flag_vector(i+1:len)=flag_vector_before(i:end);
			end
		end
	end
	if ~empty_flag
		flag_vector=flag_vector_before;
	end
end

function flag_vectors=get_flag_vectors(flag_vector_preprocessing,p);
	non_empty_flag_vector={};
	for i=1:length(flag_vector_preprocessing)
		if ~isempty(flag_vector_preprocessing{i})
			non_empty_flag_vector{end+1}=flag_vector_preprocessing{i};
		end
	end
	noMC_flag_vector={};
	for i=1:length(non_empty_flag_vector)
		if i~=p
			noMC_flag_vector{end+1}=non_empty_flag_vector{i};
		end
	end
	noMC_flag_vector_len=length(noMC_flag_vector);
	flag_vectors={};
	for i=1:(2^noMC_flag_vector_len-1)
		bi=de2bi(i);
		if length(bi)<noMC_flag_vector_len
			bi=[bi,zeros(1,noMC_flag_vector_len-length(bi))];
		end
		current_factor_vector={};
		for ii=1:length(bi)
			if bi(ii)
				current_factor_vector{end+1}=noMC_flag_vector{ii};
			else
				current_factor_vector{end+1}={};
			end
		end
		flag_vectors{end+1}=insertMC2cellArray(current_factor_vector,p,flag_vector_preprocessing);
	end
end

function num=get_NonEmpty_flag_num(flag_vector)
	num=0;
	for i=1:length(flag_vector)
		if ~isempty(flag_vector{i})
			num=num+1;
		end
	end
end

function puncta_positions=get_puncta_positions(flag_vector)
	non_empty_flag_vector={};
	for i=1:length(flag_vector)
		if ~isempty(flag_vector{i})
			non_empty_flag_vector{end+1}=flag_vector{i};
		end
	end
	puncta_positions=[];
	for i=1:length(non_empty_flag_vector)
		flag=non_empty_flag_vector{i};
		if ~isempty(flag) && strcmp(flag(1),'puncta')
			puncta_positions(end+1)=i;
		end
	end
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

function save_factor_channel_num(flag_vector)
	factor_channel_num=[];
	for c=1:length(flag_vector)
		flag=flag_vector{c};
		if ~isempty(flag) && strcmp(flag{2},'factor')
			factor_channel_num=[factor_channel_num,c];
		end
	end
	dlmwrite('factor_channel_num.txt',factor_channel_num);
end

function your_text=readTxt2cell(your_filename,numLines)
	fid = fopen(your_filename,'r');
	your_text = cell(numLines,1);
	for ii = 1:numLines
		your_text(ii) = {fgetl(fid)}; 
	end
	fclose(fid);
end