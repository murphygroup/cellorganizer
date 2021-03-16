function [f,n] = ml_removegraylevel (f,n,grayLevel)

idxGL = find(n==grayLevel);
if ~isempty(idxGL)
   f(:,idxGL)= [];
   f(idxGL,:)= [];
   n(idxGL)= [];
end