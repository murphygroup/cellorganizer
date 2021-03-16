function [uout,uerr] = unixfind( rootdir, filetype, addarg)

if ~exist('addarg','var')
    addarg = [];
end

% Perform search
[uerr uout] = unix( ['find ' rootdir ' -name "*.' filetype '" ' addarg]);

% if uerr~=0
%     error('Check inputs');
% end
