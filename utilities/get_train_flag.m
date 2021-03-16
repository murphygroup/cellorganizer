function answer = get_train_flag( options )

if isfield( options, 'train' ) && isfield( options.train, 'flag' )
    answer = options.train.flag;
end

end%get_train_flag