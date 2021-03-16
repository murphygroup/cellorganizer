function [ind1, ind2, isdone] = matrix_completion_random_columns(distances_sampled, dist_nans, options)

    options = ml_initparam(options, struct( ...
        'save_location', [], ...
        'max_columns', size(distances_sampled,2) ...
        ));
    
    save_location = options.save_location;
    max_columns = options.max_columns;

    savefile = [save_location filesep 'rand_column_recordfile.mat'];
    
    if ~isempty(dist_nans)
        distances_sampled(dist_nans) = 1;
    end

    distances_sampled = (distances_sampled+distances_sampled') > 0;
    
    if sum(sum(distances_sampled,1) == size(distances_sampled,1)) >= options.max_columns
        ind1 = -1;
        ind2 = -1;
        isdone = true;
        return
    end
    
    if exist(savefile, 'file')
        load(savefile)
    else
        working_columns_save = zeros(1, size(distances_sampled,2));
    end
    
    %find columns with more then one entry and are not entirely completed
    working_columns = working_columns_save & sum(distances_sampled,1) < size(distances_sampled,1);
    
    columninds = find(working_columns);
    
    isdone = false;
    
    %if there is a working column
    if ~isempty(columninds)
        %pick the first column that is being worked on (there should be only one in theory)
        ind2 = columninds(1);
        rowinds = find(~distances_sampled(:,ind2));
        
        %pick a random missing element in that column.
        ind1 = rowinds(randi(length(rowinds)));
        
    %if we need to select a new column
    else
        
        canstart = false;

        while ~canstart
            disp('Waiting for matrix completion file lock')
            [canstart, ~, ~, temp_name] = chunk_start([savefile  '_lock'], []);
            if ~canstart
                pause(1);
            end
        end

        working_columns = working_columns_save;
        
        if exist(savefile, 'file')
            load(savefile)
        else
            working_columns_save = zeros(size(working_columns));
        end
        
        %see if some other process initiated a new working column
        newcolumnind = find(~working_columns & working_columns_save);
        
        ind1 = 1;
        
        if ~isempty(newcolumnind)
            
            ind2 = newcolumnind(1);            
        else
            
            possible_new_columns = find(full(sum(distances_sampled,1) < size(distances_sampled,1)));
            
            if isempty(possible_new_columns)
                isdone = true;
                ind2 = -1;
            else
                new_column_ind = randi(sum(possible_new_columns>0));

                working_columns_save(possible_new_columns(new_column_ind)) = 1;

                ind2 = possible_new_columns(new_column_ind);
                
                disp('Writing updated working column file');
                save(savefile, 'working_columns_save')
            end
        end
        
        disp('Releasing working column file')
        chunk_finish(temp_name);
    end
    
    if ind1 > ind2
        tmp = ind1;
        ind1 = ind2;
        ind2 = tmp;
    end

end

