function cmtable = ml_confmat_html(cm,truesets,outsets)
%ML_CONFMAT_HTML Generate confusion matrix for HTML.
%   CMTABLE = ML_CONFMAT_HTML(CM,TRUESETS,OUTSETS) returns a string
%   showing a table of confusion matrix by HTML. CM is a confusion
%   matrix. TRUESETS is a cell array of the names of true classes and 
%   OUTSETS is a cell array of classifier outputs.

%   24-Jun-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('At least 3 arguments are required')
end

ntrueclass=length(truesets);
noutclass=length(outsets);

cmtable='<TABLE border="1">';
cmtable=[cmtable ' ' '<TR><TD rowspan=2><i>True Class</i></TD><TD align=center colspan='...
        num2str(noutclass) '><i>Output of the Classifier</i></TD></TR>'];
cmtable=[cmtable ' ' '<TR>'];
for i=1:noutclass
    cmtable=[cmtable ' ' '<TD align=center><B>' outsets{i} '</B></TD>'];
end
cmtable=[cmtable ' ' '</TR>'];

for i=1:ntrueclass
    cmtable=[cmtable ' ' '<TR>'];
    cmtable=[cmtable ' ' '<TD align=center><B>' truesets{i} '</B></TD>'];
    for j=1:noutclass
        cmtable=[cmtable ' ' '<TD align=center>' num2str(round(cm(i,j)*10)/10) '</TD>'];
    end
    cmtable=[cmtable ' ' '</TR>'];
end

cmtable=[cmtable ' ' '</TABLE>'];