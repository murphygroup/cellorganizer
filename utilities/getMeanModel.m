function [ mean_model ] = getMeanModel( models, weights )
%Compute (weighted) average model
%
%gj Aug 30, 2013

if ~exist('weights', 'var')
    weights = ones(size(models)) / length(models);
end

weights = weights ./ sum(weights);

mean_model = meanstruct(models, weights);

end

function [mean_struct] = meanstruct(struct_parents, weights)
    %Create a new struct
    mean_struct = struct;
    
    %Get all the fields
    fnames = fieldnames(struct_parents{1});
    
    %For each field
    for i = 1:length(fnames)
        %retreive the value for each field
        vals = arrayfun(@(x) eval(['struct_parents{' num2str(x) '}.' fnames{i}]), 1:length(struct_parents), 'uniformoutput', false);

        if isstruct(vals{1})
            %%%RECURSIVE%%%
            mval = meanstruct(vals, weights);
        elseif isstr(vals{1})
            %if its a string we assume all the strings are the same
            mval = vals{1};
        else
            %take the weighted average
            wval = cellfun(@(x,y) x.*y, vals, num2cell(weights), 'uniformoutput', false);
            catdim = ndims(wval{1})+1;
            wcat = cat(catdim, wval{:});
            mval = sum(wcat,catdim); 
        end
        
        eval(['mean_struct.' fnames{i} ' = mval;']);
    end 
end