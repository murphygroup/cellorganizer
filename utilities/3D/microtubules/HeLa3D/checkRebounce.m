function [flag, angle] = checkRebounce(tempmt2,imXYZ,angle_thre)

flag = true;

a = tempmt2 - imXYZ(:,end);
bs = imXYZ(:,2:end) - imXYZ(:,1:end-1);

for i = 1:size(bs,2)
    b = bs(:,i);
    angle(i) = atan2(norm(cross(a,b)),dot(a,b))/pi*180; % degree 
end

x = angle>angle_thre;
maxcon1 = max( diff( [0 (find( ~ (x > 0) ) ) numel(x) + 1] ) - 1); % find the maximum number of consecutive 1s.

if (maxcon1>=max(3,length(angle)*0.1)) || (sum(x)>=ceil(length(angle)*0.25))  %%
   flag = false;
end
