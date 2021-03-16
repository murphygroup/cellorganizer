function cropfilename = tz_getfilename_dna_2dhela(protfilename,str,ext)
%TZ_GETFILENAME_DNA_2DHELA Crop file name of 2dhela.
%   CROPFILENAME = TZ_GETFILENAME_DNA_2DHELA(PROTFILENAME) returns 
%   crop file name with key string 'dap' and extension 'dat'
%
%   CROPFILENAME = TZ_GETFILENAME_DNA_2DHELA(PROTFILENAME,STR) specifies
%   the key string by STR.
%   
%   CROPFILENAME = TZ_GETFILENAME_DNA_2DHELA(PROTFILENAME,STR,EXT)
%   specifies the key string by STR and the extension by EXT.
%   
%   See also TZ_GETFILENAME_CROP_2DHELA

%   07-Jul-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('At least 1 arguments are required')
end

if nargin<2
    str='dap';
    ext='dat';
end

if nargin<3
    ext='dat';
end

dotpos=find(protfilename=='.');
cropfilename=[protfilename(1:dotpos(1)+2) str protfilename(dotpos(2):dotpos(end)-1) '.' ext];