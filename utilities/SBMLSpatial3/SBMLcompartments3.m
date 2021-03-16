function[modelObj] = SBMLcompartments3(SBMLfile)
%This function reads an SBML file and returns the compartment list in the
%same format as sbmlimport
%
%Devin Sullivan 4/19/14

%initialize the modelObject
modelObj.Compartments = struct;


%Parse the SBML file to find compartments
docNode = xmlread(SBMLfile);
docRootNode = docNode.getDocumentElement;
allListitems = docNode.getElementsByTagName('model');
wrapperNode = allListitems.item(0);
allListitems = wrapperNode.getElementsByTagName('listOfCompartments');
items = allListitems.item(0);

%Isolate just the compartment nodes from the list 
ListOfCompartments = items.getElementsByTagName('compartment');


%This is a java struture so it uses 0 indexing!
%Loop through the compartments and add them to the structure
for i = 0:ListOfCompartments.getLength-1
    compartment = ListOfCompartments.item(i);
    compartmentID = char(compartment.getAttribute('id'));
    %use shifted index for matlab indexing
    modelObj.Compartments(i+1).name = compartmentID;
end
