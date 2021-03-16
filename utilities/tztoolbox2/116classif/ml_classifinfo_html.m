function str = ml_classifinfo_html(regmodel,option)
%ML_CLASSIFINFO_HTML Generate classifier information for html file.
%   STR = ML_CLASSIFINFO_HTML(REGMODEL) returns a string containing brief
%   summary of the trained classifier REGMODEL, which is a structure (see
%   ML_REGRESS for more details.
%   
%   STR = ML_CLASSIFINFO_HTML(REGMODEL,OPTION) returns html texts based on
%   OPTION:
%       'train' - training stage information
%       'test' - testing stage information

%   14-Jul-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('option','var')
    option = 'train';
end

switch option
case 'train'
    str=['<P>A ' ml_getmodelfullname(regmodel.modelname) ...
            ' classifier has been trained successfully.</P>'];
case 'test'
    %     if strcmp(regmodel.modelname,'lda')==1
    str = ['<P>These results are from a ' ...
        ml_getmodelfullname(regmodel.modelname) ...
        ' classifier that had an accuracy of ' ...
        num2str(round(regmodel.trconfmats.avgacc*100)) ...
        '%% on the training set'];
    if isfield(regmodel,'cvconfmats')
        str = [str ' and an accuracy of ' ...
            num2str(round(regmodel.cvconfmats.avgacc*100)) ...
            '%% using ' num2str(regmodel.cvconfmats.nfold) ...
            '-fold cross-validation'];
    end

    str = [str '.</P>'];
    
%     else
%         str=['<P> ' ml_getmodelfullname(regmodel.modelname) ...
%                 ' classifier has been trained successfully.</P>']; 
%     end
otherwise
    error('Unrecognized option.');
end

switch regmodel.modelname
case 'bpnn'
    str=[str '<P> Parameters: </P> <UL><LI>hidden nodes: ' num2str(regmodel.t.hidden) '</LI>'];
    str=[str '<LI>transfer function: logisitc</LI></UL>'];
case 'svm'
    str=[str '<P> Parameters: </P> <UL><LI>Multiclass scheme: Maxwin</LI>'];
    str=[str '<LI>kernel: rbf</LI></UL>'];
case 'lda'
    
otherwise
    error('Unrecognized classifier.');
end