function editBNGLSpatialParams3(bnglpath,sbmlfile,savepath)
%
%This function edits the bngl spatial parameter files directly using the
%sbmlspatial file. 
%This is done automatically within the sbml-spatial file, but translation
%from sbml-spatial to bngl is not yet possible.
%
%INPUTS:
%bnglpath = string containing the path to the bngl spatial parameter file
%sbmlspatialpath = string containing the path to the sbml-spatial .xml file
%savepath = string containing the path where the new bngl file should be saved
%
%OUTPUTS:
%Creates a new bngl parameter file in the location specified by savepath
%
%Created by: Devin Sullivan 11/20/14

%get the objects that need to be edited
[bnglabels,bngvalues] = textread(bnglpath,'%s%s');
tmplabels = bnglabels;
tmplabels{7} = 'sa_CP';
tmplabels{8} = 'sa_NU';
tmplabels{9} = 'sa_EN';

%get the compartment values for the existing SBML file
    docNode = xmlread(sbmlfile);
    docRootNode = docNode.getDocumentElement;
    allListitems = docNode.getElementsByTagName('model');
    wrapperNode = allListitems.item(0);
    allListitems = wrapperNode.getElementsByTagName('listOfParameters');
    ListOfParameters = allListitems.item(0);
    ParamList = ListOfParameters.getElementsByTagName('parameter');
    
    
bngvalues = cell(length(bngvalues),1);
bngvalues{1} = 'parameters';
bngvalues{end} = 'parameters';

for k = 0:ParamList.getLength-1
    thisListitem = ParamList.item(k);
    oldItem = thisListitem;
    
    % Get the label element. In this file, each
    % listitem contains only one label.
    %         thisList = thisListitem.getElementsByTagName('id');
    %         thisElement = thisList.item(0);
    idVal = thisListitem.getAttribute('id');
    
    %check to see if the parameter values correspond to one of our ids
    %fieldtoedit = strcmpi(idVal,bnglabels);
    fieldtoedit = strcmpi(idVal,tmplabels);
    %if the string does not exist, go to the next parameter field
    if sum(fieldtoedit)==0
        continue
    end
    idVal
    newparamVal = char(thisListitem.getAttribute('value'));
    
    bngvalues{fieldtoedit} = newparamVal;
    
%     % Check whether this is the label you want.
%     % The text is in the first child node.
%     if strcmp(idVal, volname)
%         thisListitem.setAttribute('value',num2str(voldata));
%         ListOfParameters.replaceChild(thisListitem,oldItem);
%         volcount = 1;
%         if sacount ==1
%             break
%         end
%         
%     elseif strcmp(idVal, saname)
%         thisListitem.setAttribute('value',num2str(sadata));
%         ListOfParameters.replaceChild(thisListitem,oldItem);
%         sacount = 1;
%         if volcount ==1
%             break
%         end
%     end
%     
    
end


%create new file 
fid = fopen(savepath,'w');
formatSpec = '%s   %s\n';
totbng = [bnglabels,bngvalues];
[nrows,ncols] = size(totbng);
for row = 1:nrows
    fprintf(fid,formatSpec,totbng{row,:});
end

%first make a table for writing out 