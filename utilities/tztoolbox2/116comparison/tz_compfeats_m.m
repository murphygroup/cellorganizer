function [pvalue,ts]=tz_compfeats_m(feats1,feats2,labeled,method,t,option)
%TZ_COMPFEATS_M Obsolete.
%
%See also ML_COMPFEATS_M

%function [pvalue,ts]=tz_compfeats_m(feats1,feats2,labeled,method,option)
%
%OVERVIEW
%   multivariate test
%PARAMETERS:
%   feats1 - feature matrix 1 n1xp
%   feats2 - feature matrix 2 n2xp
%   labeled - the feature matrices are labeled or not
%   method - multivariate test method
%   option - mixed or not
%RETURN:
%   pvalue - p-value
%   ts - test statistic
%
%HISTORY
%   18-FEB-2004 Initial write TINGZ
%   01-JUL-2004 Modified TINGZ
%       - add t for test parameters

error(tz_genmsg('of','tz_compfeats_m','ml_compfeats.m');

if ~exist('option','var')
    option='pw';
end

nansample1=find(isnan(sum(feats1,2)));
nansample2=find(isnan(sum(feats2,2)));

if ~isempty(nansample1)
    warning(['Row ' num2str(nansample1) ' in the 1st group has NaN features. Removed.']);
    feats1(nansample1,:)=[];
end


if ~isempty(nansample2)
    warning(['Row ' num2str(nansample2') ' in the 2nd group has NaN features. Removed.']);
    feats2(nansample2,:)=[];
end

if labeled==0
    comfeats=[feats1;feats2];
else
    comfeats=[feats1(:,1:(end-1));feats2(:,1:(end-1))];
end

sample_no=size(comfeats,1);
sno1=size(feats1,1);
sno2=size(feats2,1);

switch(option)
case 'pw',
    sel=1:sample_no;
case 'mx',
    sel=randperm(sample_no);
end

sel1=sel(1:sno1);
sel2=sel((sno1+1):sample_no);

s1=comfeats(sel1,:);
s2=comfeats(sel2,:);

if ~isempty(t)
    testcmd=[method '(s1,s2,' tz_cell2str(t) ')'];
else
    testcmd=[method '(s1,s2)'];
end

[pvalue,ts]=eval(testcmd);

