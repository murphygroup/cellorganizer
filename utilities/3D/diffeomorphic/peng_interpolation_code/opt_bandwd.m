function [h,H,L_h] = opt_bandwd(x)
%OPT_BANDWD choose optimal bandwidth for kernel density estimation by
%approximated logarithmic pseudo-likelihood

H = .01:.001:.5;
L_h = zeros(1,length(H));
for i = 1:length(H)
    L_h(i) = kern_est_loocv(x,H(i));
end

[maxL,maxid] = max(L_h);
h = H(maxid);
% figure
% figuresize = [8.7,6.5];
% plot(H,L_h,'k','LineWidth',2)
% xlabel('bandwidth','fontsize',6)
% ylabel('LOOCV score','fontsize',6)
% set(gca,'fontsize',fz_num)
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0,0,figuresize])
% print('-dtiff',['-r' int2str(600)],'/figure2.tif')