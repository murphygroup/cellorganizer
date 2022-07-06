function answer = check_if_SBML_output_supported( models )

answer = true;
for k=1:1:length(models)
    if isfield(models{1}, 'proteinModel') && ~strcmpi(models{1}.proteinModel.class, 'vesicle' ) && ~strcmpi(models{1}.proteinModel.type, 'gmm' )
        answer = false;
    end
end
end%check_if_SBML_output_supported
