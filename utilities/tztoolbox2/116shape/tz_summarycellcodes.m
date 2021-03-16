function s = tz_summarycellcodes(combcellcodes)
%TZ_SUMMARYCELLCODES Summarize cell codes.
%   S = TZ_SUMMARYCELLCODES(COMBCELLCODES) returns a structure containing
%   the features of the cell array of cell codes COMBCELLCODES.
%   
%   See also

%   17-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

for i=1:length(combcellcodes)
    s.cellfeats(i,:)=tz_contourfeat(combcellcodes{i}.cellcontour);
    s.arearatio(i)=combcellcodes{i}.cellarea/combcellcodes{i}.nucarea;
    s.dangle(i)=mod(combcellcodes{i}.nucmangle-combcellcodes{i}.cellmangle,180);
    
end

s.dangle(s.dangle>90)=180-s.dangle(s.dangle>90);
