function mout = findreplacestring( uout, fstr, rstr)

idx = findstr(uout, fstr);
g = [ones(1,length(fstr)) zeros(1,length(fstr))];
ind = zeros(size(uout));
ind(idx) = 1;

ind = imfilter(ind,g,'replicate');
idx2 = find(ind);

tout = uout;
tout(idx2) = [];

idx_st = 1;
mout = char(zeros(1, length(tout)+length(rstr)*length(idx)));
idx2_st = 1;
for i=1:length(idx)
    tmp_str = [uout(idx_st:idx(i)-1) rstr];
    idx2_en = idx2_st + length(tmp_str)-1;
    
    mout(idx2_st:idx2_en) = tmp_str;

    idx_st = idx(i)+length(fstr);
    idx2_st = idx2_en+1;
end
mout(idx2_st:end) = uout(idx_st:end);
