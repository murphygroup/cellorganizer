clear
load data/metadata/hpalists.mat

u = unique(classlabels)
for i=1:length(u)
    abs = unique(antibodyids(strmatch(u{i},classlabels)));
    len = length(abs);

    disp( [u{i}(5:end) ' : ' num2str(len)]);
end
