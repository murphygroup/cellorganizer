function [X,labels] = ...
    tz_transfeats_df(features,names,form,isdrug,cutoff,drugsel,protsel)
%TZ_TRANSFEATS_DF Convert DF into drug features.
%   X = TZ_TRANSFEATS_DF(FEATURES,NAMES,FORM,ISDRUG,CUTOFF,DRUGSEL,PROTSEL)
%   returns the feature matrix of drugs. FEATURES is the DF cell array with
%   NAMES for each group. FORM specifies how drug features are calculated.
%   
%   [X,LABELS] = TZ_TRANSFEATS_DF(...) aslo returns feature names.

%   16-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [X,labels] = transfeats(features,names,form,isdrug,cutoff)
%OVERVIEW:
%   convert df into drug features
%PARAMETERS:
%   features - df
%   names - group names
%   form - drug feature forms
%   isdrug - drug features or protein features
%   cutoff - threshold
%   drugsel - which drugs
%   protsel - which proteins
%RETRUN:
%   X - drug features
%   labels - feature names
%DESCRIPTION
%
%HISTORY:
%   ??-???-???? Initial write TINGZ


 
protnum=size(features,3);
drugnum=size(features,2);
featnum=size(features{1,1,1},2);
%X=zeros(drugnum,featnum*protnum);
X=[];


if ~exist('drugsel','var') | length(drugsel)==0
    drugsel=1:drugnum;
end
   

if ~exist('protsel','var') | length(protsel)==0
    protsel=1:protnum;
end
    
switch form
case 'pvalue'
    %[num_ones, pvalues, results, boths]= tz_compare_all(features,0.05);
    for drug_no=1:length(drugsel)
        for prot_no=1:length(protsel)
            both=tz_compfeats_u('tz_ttest2',features{1,drugsel(drug_no),protsel(prot_no)},features{2,drugsel(drug_no),protsel(prot_no)},0);
            if isdrug==1
                X(drug_no,(featnum*(prot_no-1)+1):featnum*prot_no)=both(2,:);
            else
                X(prot_no,(featnum*(drug_no-1)+1):featnum*drug_no)=both(2,:);
            end
        end
       % labels{drug_no}=names{2,drug_no,1}{2};
    end
case 'test1'    %Multiple test results
    for drug_no=1:length(drugsel)
        for prot_no=1:length(protsel)
            both=MultipleTestMatrix('tz_ttest2',features{1,drugsel(drug_no),protsel(prot_no)},features{2,drugsel(drug_no),protsel(prot_no)},1,0.05);
            both=sortrows(both',1)';
            if isdrug==1
                X(drug_no,(featnum*(prot_no-1)+1):featnum*prot_no)=both(2,:);
            else
                X(prot_no,(featnum*(drug_no-1)+1):featnum*drug_no)=both(2,:);
            end
        end
    end
case 'test2'    %Thresholding results
    for drug_no=1:length(drugsel)
        for prot_no=1:length(protsel)
            both=tz_compfeat('tz_ttest2',features{1,drugsel(drug_no),protsel(prot_no)},features{2,drugsel(drug_no),protsel(prot_no)},0);
            if isdrug==1
                X(drug_no,(featnum*(prot_no-1)+1):featnum*prot_no)=(both(2,:)<=0.05);
            else
                X(prot_no,(featnum*(drug_no-1)+1):featnum*drug_no)=(both(2,:)<=0.05);
            end
        end
    end
case 'diff'
    for drug_no=1:length(drugsel)
        for prot_no=1:length(protsel)
            diff=mean(features{2,drugsel(drug_no),protsel(prot_no)})-mean(features{1,drugsel(drug_no),protsel(prot_no)});
     
            if isdrug==1
                X(drug_no,(featnum*(prot_no-1)+1):featnum*prot_no)=diff;
            else
                X(prot_no,(featnum*(drug_no-1)+1):featnum*drug_no)=diff;
            end
        end
       % labels{drug_no}=names{2,drug_no,1}{2};
    end
case 'feat'
    for drug_no=1:length(drugsel)
        for prot_no=1:length(protsel)
            feat=mean(features{2,drugsel(drug_no),protsel(prot_no)});
     
            if isdrug==1
                X(drug_no,(featnum*(prot_no-1)+1):featnum*prot_no)=feat;
            else
                X(prot_no,(featnum*(drug_no-1)+1):featnum*drug_no)=feat;
            end
        end
    end    
end

if isdrug==1
    for drug_no=drugsel
        labels{drug_no}=names{2,drug_no,1}{2};
    end
else
    for prot_no=protsel
        labels{prot_no}=names{1,1,prot_no}{1};
    end
end

if cutoff<1
    i=1;
    while i<=size(X,2)
        if all(X(:,i)>cutoff)
            X(:,i)=[];
            i=i-1;
        end
        i=i+1;
    end
end