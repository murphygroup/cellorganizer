function result = as_addStack(Fi,no,val)
% AS_ADSTACK Add a 3d stack along the Z-axis, at the top and bottom
%   Fi - is the data; no is the number of stacks at the top, bottom; val is
%   the value ('1' for ones, '0' for zeros, or anything else)
%
% Created by Aabid Shariff, Murphy, Gustavo Lab

tempimg = val.*ones(size(Fi,1),size(Fi,2));

for j = 1:no
   for i = (size(Fi,3)+1):-1:2
       Fi(:,:,i) = Fi(:,:,i-1);
   end
end

for i = 1: no
   Fi(:,:,i) = tempimg;
end

for i = (size(Fi,3)+1):(size(Fi,3)+no)
   Fi(:,:,i) = tempimg;
end


result = Fi;

end