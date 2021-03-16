% This is the test for 3D Helix completion
% Author:
% Bamdev Mishra <b.mishra@ulg.ac.be>, 2011

clear all; close all; clc;

% Loading the Helix dataset
load Helix_121.mat

cd ..
% Problem formulation
trueDists = squareform(A)';
[n, ~] = size(A);
r = 3; % True dimension
p = 1; % Starting rank
f = 0.85; % Fraction of Unknown entries
params.pmax = 10; % Maximum rank


% Information
nDof = n*r - (r*(r-1)/2);
nMeas = round(f*(n*(n-1)/2));


% Compute all pair of indices
H = tril(true(n),-1);
[I,J] = ind2sub([n,n],find(H(:))); 

% Train-Test split
[train, test] = crossvalind('HoldOut',length(trueDists),f);

% Parameters
params.tol = 1e-5;
params.vtol = 1e-5;
params.verb = false;

% Initialization
Y0 = randn(n,p);


% Run algorithms

% GD
fprintf('---------- GD ----------\n');
tic;
[Y_gd infos_gd] = lowrank_dist_completion(@gd_dist_completion,I(train),J(train),trueDists(train),Y0,params);
toc;  

% TR
fprintf('---------- TR ----------\n');
tic;
[Y_tr infos_tr] = lowrank_dist_completion(@tr_dist_completion,I(train),J(train),trueDists(train),Y0,params);
toc;

% Compute predictions
estimDists_gd = pdist(Y_gd)'.^2;
estimDists_tr = pdist(Y_tr)'.^2;

% Compute mean squared errors
errorOnTest_gd = mean((trueDists(test) - estimDists_gd(test)).^2);
errorOnTrain_gd = mean((trueDists(train) - estimDists_gd(train)).^2);
errorOverall_gd = mean((trueDists - estimDists_gd).^2);

errorOnTest_tr = mean((trueDists(test) - estimDists_tr(test)).^2);
errorOnTrain_tr = mean((trueDists(train) - estimDists_tr(train)).^2);
errorOverall_tr = mean((trueDists - estimDists_tr).^2);



%% Plot results

% GD algorithm
fs = 20;
figure('name', 'Gradient descent');
line(1:length(infos_gd.costs),infos_gd.costs,'Color','r','LineWidth',2);
ax1 = gca;
set(ax1,'FontSize',fs);
xlabel(ax1,'Number of iterations','FontSize',fs);
ylabel(ax1,'Cost','FontSize',fs);
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k');
set(ax2,'FontSize',fs);
line(1:length(infos_gd.costs),infos_gd.costs,'Color','r','LineWidth',2,'Parent',ax2);               
set(ax2,'XTick',infos_gd.newRank(1:end-1),...
        'XTickLabel',2:length(infos_gd.newRank-1),...
        'YTick',[]);
    
set(ax2,'XGrid','on');
title('Rank');
legend('Gradient descent');
legend 'boxoff';




% TR algorithm
fs=20;    
figure('name', 'Trust-region');

line(1:length(infos_tr.costs),infos_tr.costs,'Marker','*','LineStyle','--','Color','blue','LineWidth',1.5);
ax1 = gca;

set(ax1,'FontSize',fs);
xlabel(ax1,'Number of iterations','FontSize',fs);
ylabel(ax1,'Cost','FontSize',fs);

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k');
       
set(ax2,'FontSize',fs);
line(1:length(infos_tr.costs),infos_tr.costs,'Marker','*','LineStyle','--','Color','blue','LineWidth',1.5,'Parent',ax2);
set(ax2,'XTick',infos_tr.newRank(1:end-1),...
        'XTickLabel',2:length(infos_tr.newRank-1),...
        'YTick',[]);
    
set(ax2,'XGrid','on');
%xlabel(ax2,'rank','FontSize',fs);
title('Rank');
legend('Trust-region');
legend 'boxoff';

% Projecting onto the 3D subspace to visualize  
rhat_gd = size(Y_gd, 2);
rhat_tr = size(Y_tr, 2);
Yhat_gd = Y_gd;
Yhat_tr = Y_tr;

if rhat_gd > 3, % we need to project onto the 3D dominant subspace
    [U1, S1, V1] = svds(Y_gd, 3);
    Yhat_gd = U1*S1*V1';
end
    
if rhat_tr > 3, % we need to project onto the 3D dominant subspace
    [U2, S2, V2] = svds(Y_tr, 3);
    Yhat_tr = U2*S2*V2';
end


% Plot 3D structure given by GD
figure('name','3D structure of GD')
fs=20;
ax1 = gca;
set(ax1,'FontSize',fs);
plot3(Yhat_gd(:,1), Yhat_gd(:,2), Yhat_gd(:,3),'*','Color', 'r','LineWidth',1.0);
xlabel(ax1,'X axis','FontSize',fs);
ylabel(ax1,'Y axis','FontSize',fs);
zlabel(ax1,'Z axis','Fontsize',fs);
legend('GD')


% Plot 3D structure given by TR
figure('name','3D structure of TR')
fs=20;
ax1 = gca;
set(ax1,'FontSize',fs);
plot3(Yhat_tr(:,1), Yhat_tr(:,2), Yhat_tr(:,3),'*','Color', 'b','LineWidth',1.0);
xlabel(ax1,'X axis','FontSize',fs);
ylabel(ax1,'Y axis','FontSize',fs);
zlabel(ax1,'Z axis','Fontsize',fs);
legend('TR')
