function filename = tz_synfilename(orf,ver,x,exposure,filter,des)
%TZ_PARSESYNFILENAME
%   FILENAME = TZ_SYNFILENAME(ORF,VER,X,EXPOSURE,FILTER) returns the
%   file name of a [yeast image].
%   
%   FILENAME = TZ_SYNFILENAME(ORF,VER,X,EXPOSURE,FILTER,DES) adds
%   coloc desginator to the file name.     
%   
%   See also tz_parseyeastfilename

%   21-Feb-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 5
    error('5 or 6 arguments are required')
end

filename = [orf '_' ver '-' x '-' exposure '_' filter];

if exist('des','var')
    if ~isempty(des)
        filename = [filename '_coloc' des];
    end
end

filename = [filename '.png'];