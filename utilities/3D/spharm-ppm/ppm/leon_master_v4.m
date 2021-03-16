[~,username] = system('whoami');
if strcmpi(strtrim(username), 'xlu2')
	root_path='/home/xlu2/'; 
	img_path_pre='/images/zhaolab/human_brain_multiplex_exm_imaging/';
elseif strcmpi(strtrim(username), 'piepi')
	root_path='/Users/piepi/';
	img_path_pre=[root_path 'CellOrganizer/Josh/'];
end
[current_path, filename, extension] = fileparts( which(mfilename) );
cd ([root_path 'CellOrganizer/cellorganizer3/']);
setup;
cd(current_path);

datetime_str=strrep(datestr(datetime),' ','_');
diary_name=[pwd filesep 'diaries' filesep 'diary_' datetime_str '.txt'];
disp(['save diary to: ' diary_name]);
diary(diary_name);
aspect=[0.161427354511474 0.161427354511474 0.4];
min_obj_size=4;
round_num={'1st' '2nd' '3rd' '4th' '5th' '6th' '7th'};

i=2; %only test for folder 2
% for i=1:length(round_num)
	disp(['running in round ' num2str(i) '...']);
	if i~=7
		imgs_paths_=[img_path_pre 'human*brain*imaging*' round_num{i} '*round/*40*.nd2' ];
	else
		imgs_paths_=[img_path_pre 'human*brain*imaging*' round_num{i} '*round/*WGA*.nd2' ];
	end
	imgs_paths=ml_ls(imgs_paths_);

	%preprocessing
	flag_vector={{'puncta','factor'},{'puncta','factor'},{'puncta','factor'},{'objects','factor'}};
	% tic;preprocessing_v3(imgs_paths,flag_vector,aspect,min_obj_size);disp('preprocessing end !');toc;
	
	%%processing all combination
	factor_vectors=get_factor_vectors();
	% p=1;
	for p=1:3
		flag_vectors={};
		% ii=7;
		for ii=1:length(factor_vectors)
			fprintf('\n%i,%i\n',p,ii);
			flag_vector=insertMC2cellArray(factor_vectors{ii},p);
			flag_vectors{end+1}=flag_vector;
			flag_vector={{'puncta','factor'},{'puncta','model'},{'puncta','factor'},{'objects','factor'}}; %this is the one cause logl>0 issue
			processing(imgs_paths,flag_vector,ii,round_num,i,root_path,img_path_pre);
			break
		end
		break
		print_best_model(flag_vectors,p,datetime_str,i)
	end

% end
print_total_run_time(diary_name);
diary off;

function processing(imgs_paths,flag_vector,ii,round_num,i,root_path,img_path_pre)
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
	tic;
	get_FactorMatrix(imgs_paths,flag_vector);disp('get factor matrix end !');
	toc;
	
	%build model
	folder_name=[img_path_pre 'human brain imaging ' round_num{i} ' round/'];
	if i~=7
		regex_input='.*40.*.nd2';
	else
		regex_input='.*WGA.*.nd2';
	end
	tic;

	[status,cmdout] = system([root_path 'anaconda/bin/python model_building_v2.py "' folder_name '" "' regex_input '" "' root_path '"'],'-echo');
	if ~isempty(strfind(cmdout,'Traceback (most recent call last):'))
		toc;
		error('Encounter error(s) when build the model !');
	else
		disp('model building end !');
		toc;
	end
end

function print_best_model(flag_vectors,c,datetime_str,i)
	diary_name=[pwd filesep 'diaries' filesep 'result_' datetime_str '_round' num2str(i) '_channel' num2str(c) '.txt'];
	disp(['save report to: ' diary_name]);
	diary(diary_name);
	full_models=dlmread('./full_model.txt');

	fid = fopen('likelihood.txt');
	likelihoods = textscan(fid,'%s%s');
	fclose(fid);
	avg=[];
	for j=1:7
		avg=[avg,str2num(likelihoods{2}{2*j-1})];
	end
	best_ind=find(avg==max(avg));

	disp('best model’s flag vector:');
	celldisp(flag_vectors(best_ind));
	disp('coefficients:');
	full_models(best_ind,:)
	disp(['likelihood:' likelihoods{2}{2*best_ind-1}]);
	disp(['likelihood sd:' likelihoods{2}{2*best_ind}]);
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

function flag_vector=insertMC2cellArray(factor_vector,p);
	len=length(factor_vector)+1;
	if p==1
		flag_vector{1}={'puncta','model'};
		flag_vector(2:len)=factor_vector;
	elseif p==len
		flag_vector{p}={'puncta','model'};
		flag_vector(1:p-1)=factor_vector;
	else
		flag_vector{p}={'puncta','model'};
		flag_vector(1:p-1)=factor_vector(1:p-1);
		flag_vector(p+1:len)=factor_vector(p:end);
	end
end

function factor_vectors=get_factor_vectors();
	factor_vectors={};
	%1 factor
	factor_vectors{end+1}={{'puncta','factor'},{},{}};
	factor_vectors{end+1}={{},{'puncta','factor'},{}};
	factor_vectors{end+1}={{},{},{'objects','factor'}};
	%2 factors
	factor_vectors{end+1}={{},{'puncta','factor'},{'objects','factor'}};
	factor_vectors{end+1}={{'puncta','factor'},{},{'objects','factor'}};
	factor_vectors{end+1}={{'puncta','factor'},{'puncta','factor'},{}};
	%3 factors
	factor_vectors{end+1}={{'puncta','factor'},{'puncta','factor'},{'objects','factor'}};
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