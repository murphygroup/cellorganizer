function idx=tz_index2order(index,msize)
%TZ_INDEX2ORDER Obsolete.
%
%See also SUB2IND

%UUU
%function tz_index2order(index,size)
%
%DESCIRPTION:
%   same as sub2ind. so no use any longer

warning('This function is obsolete and will be removed in the future');

idx=0;

index(index(:,1)>msize(1))=0;
index(index(:,2)>msize(2))=0;
index(index(:,3)>msize(3))=0;

msize=[1,cumprod(msize)];
msize=msize(1:end-1);

idx=(index-1)*msize'+1;
