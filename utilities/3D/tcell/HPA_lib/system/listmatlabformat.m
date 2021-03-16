function filelist = listmatlabformat( uout)

% Convert output to Matlab format
idx_end = findstr(uout,char(10));
idx_start = [0 idx_end(1:end-1)];

filelist = [];
filelist{length(idx_start)} = [];
for i=1:length(filelist)
    filelist{i} = uout(idx_start(i)+1:idx_end(i)-1);
end
