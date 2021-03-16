function phinm=tz_energytest2(X1,X2)
%TZ_ENERGYTEST2 Two sample multivariate energy test.
%   PHINM = TZ_ENERGYTEST2(X1,X2) returns the test statistic of the energy
%   test on two samples X1 and X2.   
%   reference:
%   A New Test for The Multivariate Two-sample Problem Based on 
%   The Concept of Minimum Energy.
%   G. Zech and B. Aslan, 2003

%   05-OCT-2003 Initial write T. Zhao
%   14-DEC-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

n=size(X1,1);
m=size(X2,1);

phia=0;
phib=0;
phiab=0;

for j=2:n
    for i=1:j-1
        phia=phia+tz_potential(X1(i,:),X1(j,:));
    end
end

phia=phia/n/n;

for j=2:m
    for i=1:j-1
        phib=phib+tz_potential(X2(i,:),X2(j,:));
    end
end


phib=phib/m/m;

for i=1:n
    for j=1:m
        phiab=phiab+tz_potential(X1(i,:),X2(j,:));
    end
end

phiab=-phiab/n/m;

phinm=phia+phib+phiab;


function dist=tz_dist(p1,p2)
dist=sqrt(sum((p1-p2).^2));

function pot=tz_potential(p1,p2)
pot=-log(tz_dist(p1,p2));
