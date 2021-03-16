function comobjfeats=tz_completeobjfeats(objfeats,classes,cellidcs,post,sel)
%TZ_COMPLETEOBJFEATS Unknown.

%function comobjfeats=tz_completeobjfeats(objfeats,classes,cellidcs,post,sel)
%
%OVERVIEW:
%   sort combined object features, cell idx and cluster label into classes and 
%PARAMETERS:
%   objfeats - combined object features
%   classes - class label
%   cellidcs - cell idx
%   post - clustering post
%   sel - feature selection
%RETURN:
%   comobjfeats - object features, cell array
%
%HISTORY:
%   24-MAR-2004 Initial write TIINGZ
%   03-NOV-2004 Modified TINGZ
%       - add comments


if size(post,2)>1
    [m,objidcs]=max(post'); 
    objidcs=objidcs';
else
    objidcs=post;
end

caclass=tz_findclass(classes);

nclass=length(caclass);

if ~isempty(sel)
    objfeats=[objfeats(:,sel),objidcs,cellidcs];
else
    objfeats=[objfeats,objidcs,cellidcs];
end
    

for i=1:nclass
    comobjfeats{i}=objfeats(classes==i,:);
end