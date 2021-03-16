function tz_objfeat_3d_p(pclass)

% img=tz_loadimage('/home/tingz/mv_3dhelaimages/ER_/cell1/prot','tif',[]);
% cropimg=imread('/home/tingz/mv_3dhelaimages/ER_/cell1/crop/ER_20020109ss12s07q1i001.mask1.tif');
% 
% objs=tz_find3dobj(img,cropimg,1);


% tz_find3dobj_mcf('/home/tingz/mv_3dhelaimages','/home/tingz/3dobjects/dnaprot','dna')
  
 features=tz_objfeat_3d_mcf('/home/tingz/3dobjectstub/prot','/home/tingz/3dobjectstub/dna','/home/tingz/3dobjects/feats5',pclass);
% features=tz_loadobjfeats_mcf('/home/tingz/3dobjects/feats');
% [combfeats,combclass,combcellidx]=tz_mcf2combobjfeats(features);
%tic
%comments=struct('time',' FEB 13, 2005','comments','results for obj paper','permID','a001');
%workdir=['/home/tingz/3dobjects/classif2']
%[avgcm,avgacc,aic]=tz_testcv(combfeats,combcellidx,combclass,[],{},10,'kmeans',19,0,'lda',1,'bpnn',0,0,'#$objnum','',workdir,workdir,workdir,[],comments,0)
%toc
