function cmtable = tz_confmat_html(cm,truesets,outsets)
%TZ_CONFMAT_HTML Obsolete

%function cmtable = tz_confmat_html(cm,truesets,outsets)
%OVERVIEW
%   
%PARAMETERS
%   cm - 
%   truesets - 
%   outsets - 
%RETURN
%   cmtable - 
%DESCRIPTION
%   
%HISTORY
%   24-Jun-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_confmat_html','ml_confmat_html'));

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