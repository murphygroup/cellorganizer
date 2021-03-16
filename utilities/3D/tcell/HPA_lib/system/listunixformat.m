function output = listunixformat( imagelist)

output = [];
for i=1:length(imagelist)
    output = [output imagelist{i} char(10)];
end

