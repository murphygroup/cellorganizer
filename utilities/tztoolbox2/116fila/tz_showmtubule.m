function tz_showmtubule(mts)
%TZ_SHOWMTUBULE Plot microtubules.
%   TZ_SHOWMTUBULE(MTS) plots microtubules in the cell array of the
%   coordinates of each filament.

%   16-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end


figure
hold on
for i=1:length(mts)
    plot(mts{i}(:,1),mts{i}(:,2));
end
hold off