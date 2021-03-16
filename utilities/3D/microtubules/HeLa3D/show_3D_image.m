%close all
%clear all

%MTpar = [150 50 0 0.995];  %% parameters of microtubules that are used for generation.
%MTpar = [150 50 0 0.95];  %% parameters of microtubules that are used for generation. %DPS 2/15/12 (.995 is not a valid colli number)
%MTpar = [150 40 0 0.95];  %% parameters of microtubules that are used for generation. %DPS 2/15/12 (50 is not a valid mulen number <=40)
MTpar = [150 15 0 0.95];  %% parameters of microtubules that are used for generation. %DPS 2/15/12 (40 is only a '-tcheck' mulen number, need real number:<=15)
%cellnum = 1;  %% the cell geometry (cell boundary and nucleus) used
cellnum = 1;  %% the cell geometry (cell boundary and nucleus) used

%[G_psf,imgcent_coordinate,imXYZ,G,mtXYZ] = getsynimage_hela(MTpar(1),MTpar(2),MTpar(3),MTpar(4),cellnum,1,1); % 
%[G_psf,imgcent_coordinate,imXYZ,G,mtXYZ] = getsynimage_hela(MTpar(1),MTpar(2),MTpar(3),MTpar(4),cellnum,1,3); %DPS 2/15/12 changed subfolder to 3 
psf = 1;%DPS 2/15/12 a boolean variable that determines whether to apply the psf to the synthesized image or not

  [G_psf,imgcent_coordinate,imXYZ,G,mtXYZ] = getsynimage_hela(MTpar(1),MTpar(2),MTpar(3),MTpar(4),cellnum,1,3,psf); %DPS 2/15/12 changed psf 
  
[protim3,Dbothfin,segdna,segcell,dnaim3,cellim3,imgcent_coordinate] = getrealimage_hela(cellnum); % 

transparency = 0.3;  %%

figure
shifted_render_data(uint8(255*double(segdna)), 'blue');
hold on
axis off
axis equal
view(-1,30)
set(gcf,'color','white');
xlim([20 220]); ylim([20 220]); zlim([20 40]);  %% to be changed for better visualization
alpha(transparency);

shifted_render_data(uint8(255*double(segcell)), [0.6,0.6,0.6]);
axis off
axis equal
view(-1,30)
set(gcf,'color','white');
xlim([20 220]); ylim([20 220]); zlim([20 40]);  %% to be changed for better visualization
alpha(transparency);

for I = 1:size(imXYZ,2), plot3(imXYZ{I}(2,:),imXYZ{I}(1,:),imXYZ{I}(3,:),'LineWidth',1.3,'Color',[1,0,0]), end
axis off
axis equal
view(-1,30)
set(gcf,'color','white');
xlim([20 220]); ylim([20 220]); zlim([20 40]);  %% to be changed for better visualization
alpha(transparency);
hold off;

saveas(gcf,'testfor.tif');