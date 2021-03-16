function answer = cleanup_options(options)
if strcmpi('nucleus',options.train.flag)
    if isfield(options,'cell')
        rmfield(options,'cell');
    end
    if isfield(options,'cell')
        rmfield(options,'protein');
    end
elseif strcmpi('cell',options.train.flag)
    if isfield(options,'cell')
        rmfield(options,'protein');
    end
end
if isfield(options,'protein')
    if !strcmpi('network',options.protein.class)
        rmfield(options.protein,'cytonuclearflag');
    end
end
answer=options;
end
