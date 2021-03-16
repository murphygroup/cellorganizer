function all_need = Find_Zeros( excel_file )
b = sum(excel_file,2);
r = find(b==0);
for i=1:length(r)-1
    if r(i+1)==r(i)+1
        r(i)=0;
    end
end
Row_=find(r~=0);
all_need=1;
for i=1:length(Row_)
    if r(Row_(i))+1 < size(excel_file, 1)
        all_need=[all_need;r(Row_(i))+1];
    end
end
end

