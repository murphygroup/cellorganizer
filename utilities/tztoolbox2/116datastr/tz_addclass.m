function all_features=tz_addclass(all_feats,feats,label)
%TZ_ADDCLASS Add a feature matrix into LMCF.
%   ALL_FEATURES = TZ_ADDCLASS(ALL_FEATS,FEATS,LABEL) addes an unlabeled 
%   feature matrix FEATS into labeled multiclass cell features(LMCF) and 
%   returns the new LMCF. LABEL is the label for the added feature matrix.
%   LABEL=-1 means a new class.

%   ??-???-???? Initial write T. Zhao
%   25-MAR-2003 Modified T. Zhao
%   13-DEC-2003 Modified T. Zhao
%   18-MAY-2004 Add comments T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

classnum=length(all_feats);
all_features=all_feats;

[samplenum featnum]=size(feats);

if(label>0)
    
    if(label>classnum) 
        warning('illegible labe'),label
        return
    end
    
    feats(:,featnum+1)=label;
    all_features{label}=[all_features{label};feats];
else
    feats(:,featnum+1)=classnum+1;
    all_features{classnum+1} = feats;
end