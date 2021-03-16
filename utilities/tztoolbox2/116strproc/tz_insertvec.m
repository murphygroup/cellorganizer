function v2=tz_insertvec(v1,iv,pos)
%TZ_INSERTVEC Inserts a vector into another vector.
%   TZ_INSERTVEC(V,IV,POS) returns a vector that is the mixture of 
%   V and IV, i.e. IV is inserted into V at the position POS. V and IV
%   must be row vectors.
%   
%   Example:
%       tz_insertvec([3 4 5],[1 2],2) returns [3 1 2 4 5]

%   19-NOV-2004 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin ~= 3
    error('Exactly 3 arguments are required');
end

if size(v1,1) > 1
    error('The 1st argument must be a row vector.');
end

if size(iv,1) > 1
    error('The 2nd argument must be a row vector.');
end

if isempty(pos)
    v2=v1;
    return;
end

spos=1;

for i=1:length(pos)
    pv2{i}=v1(spos:pos(i)-1);
    spos=pos(i);
end
pv2{i+1}=v1(spos:end);

v2=pv2{1};
for i=1:length(pos)
    v2=[v2,iv,pv2{i+1}];
end