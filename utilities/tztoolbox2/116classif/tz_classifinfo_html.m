function str = tz_classifinfo_html(regmodel)
%TZ_CLASSIFINFO_HTML Obsolete

%function str = tz_classifinfo_html(regmodel)
%OVERVIEW
%   
%PARAMETERS
%   regmodel - 
%RETURN
%   str - 
%DESCRIPTION
%   
%HISTORY
%   14-Jul-2005 Initial write TINGZ
%SEE ALSO
%

error(tz_genmsg('of','tz_classifinfo_html','ml_classifinfo_html');

str=['<P>A ' tz_getmodelfullname(regmodel.modelname) ' classifier has been trained successfully.</P>'];

switch regmodel.modelname
case 'bpnn'
    str=[str '<P> Parameters: </P> <UL><LI>hidden nodes: ' num2str(regmodel.t.hidden) '</LI>'];
    str=[str '<LI>transfer function: logisitc</LI></UL>'];
case 'svm'
    str=[str '<P> Parameters: </P> <UL><LI>Multiclass schema: Maxwin</LI>'];
    str=[str '<LI>kernel: rbf</LI></UL>'];
case 'lda'
    
otherwise
    
end