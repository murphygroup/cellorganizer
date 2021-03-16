function tz_transf_p(option)

if ~exist('hela03a_features_all','var')
  load('19990526_hela03a_all');
end

features=hela03a_features_all;
subfeats=features;
eval(['transel{1}=' option{1}]);
eval(['transel{2}=' option{2}]);

rand('state',1);
transel{3}=option{3};
[avgcm,avgacc]=tz_nfoldcvtrans_mcf(subfeats,10,'lda',{1},'pl',transel)


save(['/home/tingz/matlab/data/transf/' option{1} option{2} num2str(option{3}) '.mat']...
    ,'avgcm','avgacc');