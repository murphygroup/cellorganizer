function s = tz_strdel(s1,s2)
%TZ_STRDEL Delete substrings from another string
%   S = TZ_STRDEL(S1,S2) deletes all substrings equal to S1 in S2. S is 
%   new string.
%   Example:
%       If S1 is 'beat fate', S2 is 'e', then S is 'bat fat'

%   15-Apr-2005 Initial write TINGZ
%   Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

s=s1;
for i=1:length(s2)
    s=strrep(s,s2(i),'');
end
