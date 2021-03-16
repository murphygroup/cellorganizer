function all_features2=tz_selnamedfeats_mcf(all_features,names,selnames)
%TZ_SELNAMEDFEATS_MCF Select features by their names.
%   ALL_FEATURES2 = TZ_SELNAMEDFEATS_MCF(ALL_FEATURES,NAMES,SELNAMES)
%   returns the selected MCF according to the selected names SELNAMES.
%   NAMES is a cell array for all of the features. So the length of 
%   NAMES must be equal to the number of columns of the feature matrix.

%   ??-???-???? Initial write T. Zhao
%   02-NOV-2004 Modified T. Zhao
%       - add comments
%       - change function name tz_selnamedfeats --> tz_selnamedfeats_mcf
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

% selidx=[];
% 
% for i=1:length(selnames)
%     for j=1:length(names)
%         if strcmp(selnames(i),names(j))
%             selidx=[selidx j];
%             break;
%         end
%     end
% end

selidx = ml_strmatch(selnames,names);
for i=1:length(all_features)
    all_features2{i}=double(all_features{i}(:,selidx));
end
