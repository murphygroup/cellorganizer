function machine = tz_machineinfo(varargin)
%TZ_MACHINEINFO Get information of the machine
%   TZ_MACHINEINFO(INFO1,...,INFON) returns machine information, which is 
%   stored in a struecture.
%
%   Options for information:
%       'name': machine name
%       'computer': computer type (for details, type 'help computer')

%   24-May-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

for i=1:length(varargin)
    switch varargin{i}
        case 'name'
            if exist('machine.mat','file')
                load machine.mat
            else
                [s,hostname]=unix('hostname');
                if uint8(hostname(end))==10
                    hostname=hostname(1:end-1);
                end
                dotpos=find(hostname=='.');
                if isempty(dotpos)
                    machine.name=hostname;
                else
                    machine.name=hostname(1:dotpos-1);
                end
            end
        case 'computer'
            machine.computertype = computer;
        case 'domain'
            [s,name]=unix('dnsdomainname');
            if uint8(name(end))==10
                name = name(1:end-1);
            end
            machine.domainname = name;
    end
end