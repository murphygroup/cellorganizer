function Ixy = mutualinformation(I, NBINS)

tol = 1e-3;

s = size(I);
if s(2)~=2,
    error('improper matrix size');
end

mi = double(min(I,[],1));
ma = double(max(I,[],1));
% I = double(I);
for i=1:s(2)
    I(:,i) = (I(:,i)-mi(i)) * ((NBINS(i)-1)/(ma(i)-mi(i)));  % this is fast
    % I(:,i) = round((I(:,i)-mi(i))/(ma(i)-mi(i))*(NBINS(i)-1));  % this is exact
end

ma = double(max(I,[],1))+1;
rma = cumprod([1 ma]);

cjoint = zeros(ma);

idx = ones(s(1),1);
for i=1:s(2)
    idx = idx + double(I(:,i))*rma(i);
end

for i=1:s(1)
    cjoint(idx(i)) = cjoint(idx(i))+1;
end

cx = sum(cjoint,2);
cy = sum(cjoint,1);

cx = max(cx,tol);
cy = max(cy,tol);
cjoint = max(cjoint,tol);

Ixy = cjoint.*log(s(1)*cjoint./repmat(cx,[1 ma(2)])./repmat(cy,[ma(1) 1]));
Ixy = sum(Ixy(:))/s(1);

return
