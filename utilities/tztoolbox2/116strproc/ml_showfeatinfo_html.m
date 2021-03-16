function htmlinfo = ml_showfeatinfo_html(setname,setidx,varargin)
%ML_SHOWFEATINFO_HTML Show feature information of the features.
%   HTMLINFO = ML_SHOWFEATINFO_HTML(SETNAME,SETIDX) returns a string that
%   represents a talbe of feature information for HTML file. The features
%   are specified by a string SETNAME and an integer vector SETIDX.
%   
%   See also

%   20-Dec-2005 Modified from Sam's code  T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('At least 2 arguments are required')
end

[featdesp,featsetnames,featids] = ml_featinfo;

featureSetNameIdx = strmatch(setname,featsetnames,'exact');

selectedFeatureIds = featids{featureSetNameIdx}(setidx);


htmlinfo = '<HTML>\n';
htmlifno = [htmlinfo '<TABLE align=center cellpadding=2><TBODY>\n'];

htmlinfo = [htmlinfo '<TABLE align=center cellpadding=2><TBODY>\n'];
htmlinfo = [htmlinfo '<TD><B><U>Feature ID</B></U></TD><TD><B><U>Feature Name</B></U></TD><TD><B><U>Feature Description</B></U></TD>'];
    
if nargin>2
    colname = varargin{1};
    colvalues = varargin{2};
    htmlinfo = [htmlinfo '<TD><B><U>' colname '</B></U></TD>'];
end

htmlinfo = [htmlinfo '\n'];

for i=1:length(setidx)
    htmlinfo =  ...
        [htmlinfo '<TR><TD>' num2str(selectedFeatureIds(i)) '</TD>'];
    featureId = selectedFeatureIds(i);
    featureDescription = featdesp{featureId};
    featurelink = ml_featlink(featureId,featureDescription);
    slfName = [setname '.' num2str(setidx(i))];
    htmlinfo = ...
        [htmlinfo '<TD><a href=' featurelink ' target="_blank">' slfName '</a></TD>'];
    htmlinfo = [htmlinfo '<TD>' featureDescription '</TD>'];
    if nargin>2
       htmlinfo = [htmlinfo '<TD>' num2str(colvalues(i)) '</TD></TR>'];
    end
    htmlinfo = [htmlinfo '\n'];
end

htmlinfo =  [htmlinfo '</TBODY></TABLE></BODY>\n'];

%%
htmlinfo = [htmlinfo '</html>'];
