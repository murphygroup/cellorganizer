function [avgcm,avgacc,aic]=tz_execfold(fold)

%features=tz_loadobjfeats_mcf('/home/tingz/3dobjects/feats5');
%[combfeats,combclass,combcellidx]=tz_mcf2combobjfeats(features);
load /home/tingz/3dobjects/combfeats.mat

workdir=['/home/tingz/3dobjects/classif4']
%keyboard
[avgcm,avgacc,aic]=tz_objclassif(combfeats,combcellidx,combclass,[],{},[],[],...
    10,fold,'kmeans',19,[],'dist',{},'bpnn',0,{},'#!objnum',workdir,6,...
    {[],'leavetrain'})
