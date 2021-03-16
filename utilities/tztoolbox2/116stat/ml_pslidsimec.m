function ml_pslidsimec(featfile1,featfile2,conf,testmeth,resultfile, ...
    setname1,setname2,featuresetname,host,testname)
%ML_PSLIDSIMEC Do comparison for PSLID.
%   ML_PSLIDSIMEC(FEATFILE1,FEATFILE2,CONF,TESTMETH,RESULTFILE,SETNAME1,
%       SETNAME2,FEATURESETNAME,HOST,TESTNAME) is a function for PSLID
%       SImEC service. Here are the descriptions of the arguments:
%   FEATFILE1 - feature file path of the first set
%   FEATFILE2 - feature file path of the second set
%   CONF - confidence level
%   TESTMETH - testing method
%   RESULTFILE - file path of result report
%   SETNAME1 - name of the first set
%   SETNAME2 - name of the second set
%   FEATURESETNAME - name of the feature set
%   HOST - website for supplementary description
%   TESTNAME - standard name of the testing method
%   
%   See also

%   07-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 10
    error('Exactly 10 arguments are required')
end

[features1,featids,slfnames,featnames,featdesps] = ...
    ml_loadbinfeat(featfile1);
[features2,featids2] = ...
    ml_loadbinfeat(featfile2);

if isempty(featids)
    featids = featids2;
end

conf = round(conf*100)/100;
succ=1;
pvalue=1;
[pvalue,both]=ml_multest2(double(features1),double(features2), ...
    'Un','ttest2',1,conf,[]);
[pvalue,ts]=ml_compfeats_m(double(features1),double(features2),testmeth);
ml_simec_report(resultfile,pvalue,both,features1,features2, ...
    setname1,setname2,featuresetname,conf,succ,host,testname,featids);
