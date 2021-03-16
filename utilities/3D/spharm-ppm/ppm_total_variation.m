function [ vol_norm ] = ppm_total_variation( models )
%Total variation given two models
model1=models{1};
model2=models{2};

%MODEL1
a1=strsplit(model1.proteinModel.ppm.full.param);
 m1=str2double(a1(1:size(a1,2)-1));

 a1_min=strsplit(model1.proteinModel.ppm.full.covars_min);
 m1_min=str2double(a1_min(1:size(a1_min,2)-1));

 a1_max=strsplit(model1.proteinModel.ppm.full.covars_max);
 m1_max=str2double(a1_max(1:size(a1_max,2)-1));

 %MODEL2
 a2=strsplit(model2.proteinModel.ppm.full.param);
 m2=str2double(a2(1:size(a2,2)-1));

 a2_min=strsplit(model2.proteinModel.ppm.full.covars_min);
 m2_min=str2double(a2_min(1:size(a2_min,2)-1));

 a2_max=strsplit(model2.proteinModel.ppm.full.covars_max);
 m2_max=str2double(a2_max(1:size(a2_max,2)-1));

 %MIN 
 mn=zeros(2,1);
 mx=zeros(2,1);
 for i=1:2
  mn(i)=min(m1_min(i),m2_min(i));
  mx(i)=max(m1_max(i),m2_max(i));
 end

 if(m1==m2)
   vol_norm=1;
 else
   xdata = mn(1)+[0; rand(100,1); 1]*(mx(1)-mn(1));
   ydata = mn(2)+[0; rand(100,1); 1].*(mx(2)-mn(2));
   x = sort(xdata);
   y = sort(ydata);

  [X,Y] = meshgrid(x,y);
  Z1 = m1(1)+ m1(2)*X+ m1(3)*Y;
  vol_1=trapz(y,trapz(x,Z1,2),1);
  Z1_norm=Z1./vol_1;

  Z2=  m2(1) + m2(2)*X + m2(3)*Y;
  vol_2=trapz(y,trapz(x,Z2,2),1);
  Z2_norm=Z2./vol_2;

  Zshared=min(Z1_norm,Z2_norm);

  vol_norm=trapz(y,trapz(x,Zshared,2),1);
 end


end

