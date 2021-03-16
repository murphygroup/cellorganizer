function caclass=tz_findclass(class,ispint)
%TZ_FINDCLASS Obsolete.
%
%See also ML_FINDCLASS

%function caclass=tz_findclass(class,ispint)
%OVERVIEW
%   reorganize classes
%PARAMETERS:
%   class - class labels or labeled data with a class label in the last column
%   ispint - whether labels are positive integers
%       default: positive integer
%RETRUN:
%   class - a cell array of classes
%DESCRIPTION:
%
%HISTORY
%   APR-18-2004 Initial write TINGZ

error(tz_genmsg('of','tz_findclass','ml_findclass'));

if ~exist('ispint','var')
    ispint=1;
end

if ispint==0
  
    caclass={};
    caclass{1}=[class(1,:),1];
    nsample=size(class,1);
    assigned=0;
    
    for i=2:nsample
        assigned=0;
        for j=1:length(caclass)
            if (class(i,end)==caclass{j}(1,end-1))
                caclass{j}=[caclass{j};[class(i,:),i]];
                assigned=1;
            end
        end
        
        if assigned==0
            caclass{length(caclass)+1}=[class(i,:),i];
        end
    end
else
    label=class(:,end);
    post=ml_label2post(label);
    nclass=max(label);
    if any(sum(post,1)==0)
        post=tz_stdclst(post);
        label=ml_post2label(post);
        nclass=max(label);
    end
    samples=1:size(class,1);
    
    %class label and sample indices
    for i=1:nclass
        caclass{i}=[class(label==i,:),samples(label==i)'];
    end
end