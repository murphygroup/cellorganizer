%%%%%%%%%%%%%%%% HELPER METHOD - ADD DOCUMENTATION TO MODEL %%%%%%%%%%%%%%%
function model = add_documentation_to_model( model, param )
try
    model.documentation = param.documentation;
    model.documentation.date = date;
catch
    model.documentation.date = date;
end
end%add_documentation_to_model
