function setlist = ml_hashfeatset(featureid)
%ML_HASHFEATSET Get feature set names from feature ids.
%   SELLIST = ML_HASHFEATSET(FEATUREID) returns a cell array of strings
%   representing the found feature set names. Only feature names containing
%   any id in the vector FEATUREID will be included. If there is any ID
%   in FEATUREID is not recognized, a string 'unknown' will be included.
%   Currently the available feature sets are
%
%   {'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg' 'wav' 'gab'}
%   
%   See ML_FEATURES for more details.


%   19-Jul-2005 Initial write TINGZ T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

allfeatset={'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg' 'wav' 'gab'};
allfeattable=[ones(5,1);ones(1,1)+1;ones(14,1)+2;ones(3,1)+3; ...
        ones(49,1)+4;ones(13,1)+5;ones(5,1)+6;ones(30,1)+7;ones(60,1)+8];

unknown=0;
if any(featureid>length(allfeattable))
    unknown=1;
    featureid(featureid>length(allfeattable))=[];
end
if any(featureid<=0)
    unknown=1;
    featureid(featureid<=0)=[];
end

featsetindex=allfeattable(featureid);
flags=zeros(length(allfeatset),1);

flags(featsetindex)=1;

setlist=allfeatset(find(flags==1));

if unknown==1
    setlist={setlist{:} 'unknown'};
end