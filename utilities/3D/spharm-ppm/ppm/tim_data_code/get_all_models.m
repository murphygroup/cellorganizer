function answer = get_all_models(imgs_paths,masks_paths,intermediate_output_path,final_output_path,python_path,flag_vector,options)
% get_all_models generates all models with all possible factor combinations of all images. The resulting models and log-likelihoods will be saved in final_output_path. For example, if the flag_vector argument is set as {{},{'object',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}}}, the function will generate six models
% 1. A model of puncta 1 (channel 3) using one factor: distance to object found using channel 2;
% 2. A model of puncta 1 (channel 3) using one factor: distance to puncta found using channel 4;
% 3. A model of puncta 1 (channel 3) using two factors: distance to object found using channel 2 and distance to puncta found using channel 4;
% 4. A model of puncta 2 (channel 4) using one factor: distance to object found using channel 2;
% 5. A model of puncta 2 (channel 4) using one factor: distance to puncta found using channel 3;
% 6. A model of puncta 2 (channel 4) using two factors: distance to object found using channel 2 and distance to puncta found using channel 3.
%
% Leave-one-out cross-validation will be done in cell level when image masks are given. If image masks are not given, cross-validation will be done in image level.
%
% Input
% -----
% imgs_paths: A cell of strings that constains images' paths.
% masks_paths: A cell of strings that constains masks' paths. The order in masks_paths should match the order in imgs_paths. It can be set as empty, then the cross-validation will be done in image level.
% intermediate_output_path: The path that contains distance transformation matrices resulted from function get_all_distance_transformation_matrix. It is used to save the intermediate files. This argument should be the same intermediate_output_path as in function get_all_distance_transformation_matrix.
% final_output_path: The path that saves the resulting log-likelihood and model text file.
% python_path: The python binary path. It will be used for running model building module.
% flag_vector: A three dimension cell containing flags. It can contains
%				0 or 1 from
%					puncta = find puncta for this channel
%					objects = find objects for this channel
%					multiple_objects = find multiple objects for this channel
%					intensity = remain image intensity for this channel
%				0, 1 or 2 from
%					model = build model for the centers after preprocessing this channel
%					factor = use as factor after preprocessing and distance transform
%				e.g., {{},{'object',{'factor'}},{'puncta',{'model','factor'}},{'puncta',{'model','factor'}}} means
%				1. build a model of puncta 1 (channel 3) using one factor: distance to object found using channel 2;
%				2. build a model of puncta 1 (channel 3) using one factor: distance to puncta found using channel 4;
%				3. build a model of puncta 1 (channel 3) using two factors: distance to object found using channel 2 and distance to puncta found using channel 4;
%				4. build a model of puncta 2 (channel 4) using one factor: distance to object found using channel 2;
%				5. build a model of puncta 2 (channel 4) using one factor: distance to puncta found using channel 3;
%				6. build a model of puncta 2 (channel 4) using two factors: distance to object found using channel 2 and distance to puncta found using channel 3.
%				Note: only puncta channel can be set as model channel.
%				Channels that have empty flags mean they will be ignored.
%				This argument must be the same flag_vector as in function
%				get_all_distance_transformation_matrix. % Shen Jin Sep 04,
%				2019: add a new type "objects_local", which will use local
%				threshold to extract objects and replace it using the
%				centers
% options.dummy_num: Number of dummny points used in build_model.
% options.rand_num: Number of random points used in build_model.
% options.cv_mode: Two cross-validation modes, one is rd_img (random fold assignment in image level), the other is rd_roi (random fold assignment in ROI level).
% options.fold: Number of fold used for cross-validation in build_model.
% options.cv_round: Number of round for cross-validation in build_model.
% options.debug: A boolean value indicates whether visualize the factor matrix. The resulting figures will save in intermediate_output_path/debug folder.
% options.mask_inverted_color_flag: A boolean value indicates whether mask color need to invert.
%
% Output
% ------
% answer: A boolean indicates whether the function successfully finished or not.

% Author: Xin Lu
%
% Copyright (C) 2012-2019 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
answer=1;
% try
%loop each model channel
model_channels=get_model_channels(flag_vector);
for mc_ind=1:length(model_channels)
    mc=model_channels{mc_ind};
    % 		tic;get_FactorMatrix(imgs_paths,masks_paths,intermediate_output_path,flag_vector,mc,options);disp('get factor matrix end !');toc;
    
    tic;
    disp('Getting factor matrix');
    get_FactorMatrix_auas(imgs_paths,masks_paths, ...
        intermediate_output_path,flag_vector,mc,options);
    toc;
    
    % processing all possible flags
    flag_vectors=get_flag_vectors(flag_vector,mc);
    for ii=1:length(flag_vectors)
        % ii=length(flag_vectors); %%% need to comment out
        fprintf('\n%i,%i\n',mc,ii);
        build_model(intermediate_output_path,final_output_path, ...
            python_path,flag_vectors{ii},ii,options,mc);
        % break
    end
    % break
    print_best_model(flag_vectors,mc,options.datetime_str,final_output_path,imgs_paths,length(flag_vectors));
    % break
end
% catch
% 	answer=0;
% 	warning('get_all_models has error');
% end
end

function build_model(intermediate_output_path,final_output_path,python_path,flag_vector,ii,options,mc)
%Get factor matrix and build model

%celldisp(flag_vector);
if ~exist(final_output_path)
    system(['mkdir -p ' final_output_path]);
end
if ii==1
    if exist([final_output_path filesep 'full_model.txt'])
        delete([final_output_path filesep 'full_model.txt']);
    end
    if exist([final_output_path filesep 'likelihood.txt'])
        delete([final_output_path filesep 'likelihood.txt']);
    end
end
save_factor_channel_num(flag_vector);

%build model
tic;
command = [python_path ' ' which('model_building.py') ' "' final_output_path '" "' intermediate_output_path '" "' options.dummy_num '" "' options.rand_num '" "' options.cv_mode '" "' options.fold '" "' options.cv_round '" "' int2str(mc) '"'];
[status,cmdout] = system( command,'-echo' );
if ~isempty(strfind(cmdout,'Traceback (most recent call last):'))
    toc;
    error('Encounter error(s) when build the model!');
else
    disp('Finished buildig model');
    toc;
end
end

function model_channels=get_model_channels(flag_vector)
model_channels={};
for i=1:length(flag_vector)
    flag=flag_vector{i};
    if ~isempty(flag)
        if any(strcmp(flag{2},'model'))
            model_channels{end+1}=i;
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

function flag_vectors=get_flag_vectors(flag_vector,mc);
%flag_vector to flag_vector_preprocessing(in old version)
flag_vector_preprocessing={};
for i=1:length(flag_vector)
    flag=flag_vector{i};
    if isempty(flag)
        flag_vector_preprocessing{end+1}={};
    elseif i==mc
        flag_vector_preprocessing{end+1}={flag{1},'model'};
    elseif any(strcmp(flag{2},'factor'))
        flag_vector_preprocessing{end+1}={flag{1},'factor'};
    else
        flag_vector_preprocessing{end+1}={};
    end
end
% celldisp(flag_vector_preprocessing);
%mc to p
p=mc;
for i=1:length(flag_vector_preprocessing)
    if i==mc
        break
    elseif isempty(flag_vector_preprocessing{i})
        p=p-1;
    end
end

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
    bi=dec2bin(i);
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
flag_vector=flag_vector_preprocessing;
ii=0;
for i=1:length(flag_vector)
    if ~isempty(flag_vector{i})
        ii=ii+1;
        if isempty(flag_vector_before{ii})
            flag_vector{i}={};
        end
    end
    
end
end

function print_best_model(flag_vectors,mc,datetime_str,final_output_path,imgs_paths,flag_vector_length)
[img_folder_name,name,ext]=fileparts(imgs_paths{1});
img_folder_name=strrep(img_folder_name,filesep,'-');
% 	diary_name=[final_output_path filesep 'result_' datetime_str '_folder-' img_folder_name '_ModelChannel-' num2str(mc) '.txt'];
diary_name=[final_output_path filesep 'result_' datetime_str '_ModelChannel-' num2str(mc) '.txt'];
disp(['save report to: ' diary_name]);
diary(diary_name);

numLines=length(flag_vectors);
txt_filename=[final_output_path filesep 'full_model.txt'];
full_models=readTxt2cell(txt_filename,numLines);

fid = fopen([final_output_path filesep 'likelihood.txt']);
likelihoods = textscan(fid,'%s%s');
fclose(fid);
avg=[];
for j=1:flag_vector_length
    avg=[avg,str2num(likelihoods{2}{2*j-1})];
end
best_ind=find(avg==max(avg));

disp('Best model vector:');
celldisp(flag_vectors(best_ind));
disp(['Coefficients: ' full_models{best_ind}]);
disp(['Likelihood: ' likelihoods{2}{2*best_ind-1}]);
disp(['Likelihood SD: ' likelihoods{2}{2*best_ind}]);

diary off;
end

function your_text=readTxt2cell(your_filename,numLines)
fid = fopen(your_filename,'r');
your_text = cell(numLines,1);
for ii = 1:numLines
    your_text(ii) = {fgetl(fid)};
end
fclose(fid);
end