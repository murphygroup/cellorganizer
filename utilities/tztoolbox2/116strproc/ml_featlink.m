function link = ml_featlink(slfid,slfdesp)
%ML_FEATLINK Web link of the features.
%   LINK = ML_FEATLINK(SLFID,SLFDESP) returns a string that is 
%   the web link of the feature that has ID SLFID. If SLFID is
%   string, it is a feature set name.
%
%   See also

%   20-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

webserver = 'http://murphylab.web.cmu.edu/services/SLF/'; 

if isstr(slfid)
   featlink = [];
else
   featlink = ['features.html' '#' num2str(slfid)];

   if slfid>23.5 & slfid<72.5
      idlink='Z';
      featlink = ['Z_ims/' slfdesp '.html'];
   end
end

link = [webserver featlink];
