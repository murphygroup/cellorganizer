function Ixy = mutualinformation(I)

tol = 1e-3;

s = size(I);
if s(2)~=2,
    error('improper matrix size');
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
