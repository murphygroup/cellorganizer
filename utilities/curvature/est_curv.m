function curvature = est_curv(obj,co,i,N,KNN_half)

index = i-KNN_half:i+KNN_half;

%if it is at begin
ind = index<1;
index(ind) = index(ind)+N;
%if it is at end
ind = index>N;
index(ind) = index(ind)-N;

data = co(index,:);

[p,S,mu] = polyfit(data(:,1),data(:,2),1);

[z, r, residual] = fitcircle(data');



if S.normr > residual
        curvature = 1/r;
    
else
    curvature = 0;
end

