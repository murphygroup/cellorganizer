function setlist = tz_hashfeatset(featureid)
%TZ_HASHFEATSET Obsolete
%
%See also ML_HASHFEATSET

%function setlist = tz_hashfeatset(featureid)
%OVERVIEW
%   get set list from feature ids
%PARAMETERS
%   featureid - a vector of feature ids
%RETURN
%   setlist - list of sets
%DESCRIPTION
%   
%HISTORY
%   19-Jul-2005 Initial write TINGZ
%SEE ALSO
%   

error('of','tz_hashfeatset','ml_hashfeatset');

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