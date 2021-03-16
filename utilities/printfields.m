function printfields (structure, level, prefx)

pretab = [];
if exist('level','var')
    for i=1:level
        pretab = [pretab '>'];
    end
%    disp(pretab)
else
    level = 0;
    prefx = [];
end

try
    siz=size(structure);
    if prod(siz)>50
        whos structure
    else
        disp(structure)
    end
catch
    disp(structure)
end
try
    fields = fieldnames(structure);
    for i=1:length(fields)
%        fprintf('%s %s %s\n',prefx,fields{i},pretab)
        prefxi=[prefx '.' fields{i}];
        eval(['printfields(structure.' fields{i} ',' num2str(level+1) ...
            ', prefxi )'])
    end
catch
end
