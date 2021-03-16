function [avgcm,avgacc,ncm,kappa,ua,pa]=tz_calcavgcm(cvcm,merge)
%TZ_CALCAVGCM Obsolete.

%function [avgcm,avgacc]=tz_calcavgcm(cvcm)
%
%OVERVIEW:
%   Calculate average confusion matrix and average accuracy from cross validation results
%PARAMETERS:
%   cvcm - confusion matrices from cross validatioin; cell of matrices with numbers
%   merge - indices of classes to be merged
%RETRUN:
%   avgcm - average confusion matrix; percentage
%   avgacc - average accuracy
%DESCRIPTION:
%   Calculate confusion matrix by averaging numbers
%HISTORY:
%   17-MAY-2004 Modified TINGZ
%       - change parameter meaning

error(tz_genmsg('of','tz_calcavgcm','ml_calcavgcm'));

if ~exist('merge','var')
    merge=[];
end

%initialization
for i=1:length(cvcm)
    if ~isempty(cvcm{i})
        ncm=zeros(size(cvcm{i}));
        break;
    end
end

%add confusion matrices of all folds
for i=1:length(cvcm)
    if ~isempty(cvcm{i})
        ncm=ncm+cvcm{i};
    end
end

%merge classes
if ~isempty(merge)
    ncm2=ncm;
    ncm2(merge(1),:)=sum(ncm(merge,:),1);
    ncm2(merge(2:end),:)=[];
    
    ncm2(:,merge(1))=sum(ncm2(:,merge),2);
    ncm2(:,merge(2:end))=[];
    
    ncm=ncm2;
end

%calculate accurate rate
avgcm=ncm./repmat(sum(ncm,2),1,size(ncm,2))*100;
avgacc=sum(diag(ncm))/sum(ncm(:));

srow=sum(ncm,2);
scol=sum(ncm,1)';

%users accuracy and producers accuracy
ua=(diag(ncm)./scol)';
pa=(diag(ncm)./srow)';

%kappa
n=sum(ncm(:));
kappa=(n*sum(diag(ncm))-sum(srow.*scol))/(n^2-sum(srow.*scol));