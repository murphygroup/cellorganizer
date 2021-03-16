function [ cellparam ] = model2param( model, options )
%MODEL2PARAM Summary of this function goes here
%   Detailed explanation goes here

%topologically sort the dependency graph
synthorder = graphtopoorder(model.dependencies);

%for each component from top down
for i = 1:length(synthorder)
    %the components parameters come first, then the parameters as follows
    %in numerical order in param.dependencies
    
    %get the name of the component
    component_num = synthorder(i);
    component_name = model.component_names(component_num);
    
    
    parents = find(model.dependencies(component_num,:));
    
    parent_parameters = cell(1,length(parents));
    for j = 1:length(parents)
        parent_parameters{model.dependencies(component_num,parents(j))} = cellparam.(parents(j));
    end
    
    cellparam.(component_name) = draw_sample(model.(component_name), parent_parameters{:});
    
end


end

