function [vol_norm]= total_variation_joint_pca(X1,X2,coeff)

score1=X1*coeff;
score2=X2*coeff;
vol_norm=0;
if isequal(score1,score2)
    vol_norm=1;
else 
 [m1,s1]=normfit(score1(:,1:2));
 [m2,s2]=normfit(score2(:,1:2));
 x1=mvnrnd(m1,s1,100);
 x2=mvnrnd(m2,s2,100);
 y1=mvnpdf(x1,m1);
 y2=mvnpdf(x2,m2);
 xi = linspace(min(min(x1(:,1),x2(:,1))),max(max(x1(:,1),x2(:,1))),100);
 yi = linspace(min(min(x1(:,2),x2(:,2))),max(max(x1(:,2),x2(:,2))),100);
 dx = xi(2) - xi(1);
 dy = yi(2) - yi(1);


 [X, Y] = meshgrid(xi,yi);
 Z1=griddata(x1(:,1),x1(:,2),y1,X,Y);
 dv = Z1(~isnan(Z1))*dx*dy;
 vol_1=sum(dv);
 Z1_norm=Z1/vol_1;

 Z2=griddata(x2(:,1),x2(:,2),y2,X,Y);
 dv = Z2(~isnan(Z2))*dx*dy;
 vol_2=sum(dv);
 Z2_norm=Z2/vol_2;

 Z_shared=min(Z1_norm,Z2_norm);
 dv = Z_shared(~isnan(Z_shared))*dx*dy;

 vol_norm=sum(dv);
end
end