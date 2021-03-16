function [docNode,wrapperNode] = addVol_SurfArea3(docNode,wrapperNode,objname,voldata,sadata)

%Grab the parameters
allListitems = wrapperNode.getElementsByTagName('listOfParameters');
ListOfParameters = allListitems.item(0);
if isempty(ListOfParameters)
    ListOfParameters = docNode.createElement('listOfParameters');
    paramval = docNode.createElement('parameter');
    paramval.setAttribute('id','PlaceHolder');
    paramval.setAttribute('value','PlaceHolder');
    
end
ParamList = ListOfParameters.getElementsByTagName('parameter');

%     %find the old volume and surface area nodes
%     oldVol = ListOfParameters.getChild('vol');
%     oldSA = ListOfParameters.getChild('sa');
%
volcount = 0;
sacount = 0;
volname = ['vol_',objname];
saname = ['sa_',objname];
for k = 0:ParamList.getLength-1
    thisListitem = ParamList.item(k);
    oldItem = thisListitem;
    
    % Get the label element. In this file, each
    % listitem contains only one label.
    %         thisList = thisListitem.getElementsByTagName('id');
    %         thisElement = thisList.item(0);
    idVal = thisListitem.getAttribute('id');
    
    
    % Check whether this is the label you want.
    % The text is in the first child node.
    if strcmp(idVal, volname)
        thisListitem.setAttribute('value',num2str(voldata));
        ListOfParameters.replaceChild(thisListitem,oldItem);
        volcount = 1;
        if sacount ==1
            break
        end
        
    elseif strcmp(idVal, saname)
        thisListitem.setAttribute('value',num2str(sadata));
        ListOfParameters.replaceChild(thisListitem,oldItem);
        sacount = 1;
        if volcount ==1
            break
        end
    end
    
    
end

if volcount==0
    %make the volume node
    disp(['Volume parameter not found, appending parameter named: ',volname]);
    volume = docNode.createElement('parameter');
    volume.setAttribute('id',['vol_',objname]);
    volume.setAttribute('value',num2str(voldata));
    volume.setAttribute('constant','true');
    ListOfParameters.appendChild(volume);
end

if sacount==0
    %make the surface area node
    surface_area = docNode.createElement('parameter');
    surface_area.setAttribute('id',[saname]);
    surface_area.setAttribute('value',num2str(sadata));
    surface_area.setAttribute('constant','true');
    disp(['Surface area parameter not found, appending parameter named: ','sa_',objname]);
    ListOfParameters.appendChild(surface_area);
end

%append the data back on the wrapper
wrapperNode.appendChild(ListOfParameters);