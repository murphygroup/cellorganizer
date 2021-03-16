function ts = tz_strtok_pos(s,pos)
%TZ_STRTOK_POS Extract tokens in string.
%   TS = TZ_STRTOK_POS(S,POS) returns all tokens in the string 
%   S according to positions POS. POS is an increasing integeter 
%   vector or string. The interger vector specifies position. and
%   the string specifies token patterns.
%   TS is a cell array of these strings.
%
%   Example:
%       tz_strtok_pos('Enjoyusingthistoolbox!',[5 10 14]) returns 
%       {'Enjoy','using','this','toolbox!'}
%
%       tz_strtok_pos('Enjoyusingthistoolbox!','11111000001111') 
%       returns {'Enjoy','using','this','toolbox!'}
%
%   See also TZ_STRTOK

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if ischar(pos)
    pos = pos(2:end) - pos(1:end-1);
    if length(pos) < length(s)-1
        pos(end+1) = 1;
    end
    pos = find(pos);
end

if isempty(pos)
    ts{1} = s;
    return;
end

ts = {};

if (pos(1) >= 1)
    ts{1} = s( 1 : pos(1) );
end

if length(pos) > 1
    for i = 1:length(pos)-1
        if (pos(i) >= 1) & (pos(i) < length(s) )
            if pos(i+1) > length(s)
                pos(i+1) = length(s);
            end
            ts{end+1} = s( pos(i)+1:pos(i+1) ); 
        end
    end
end

if ( pos(end) < length(s) )
    ts{end+1} = s(pos(end)+1:end);
end

