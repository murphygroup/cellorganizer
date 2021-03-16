function listOfSpecies = addSpecies(docNode)

listOfSpecies = docNode.createElement('listOfSpecies');
species = docNode.createElement('species');
% cellorg = docNode.createElement('cellorganizer');
listOfSpecies.appendChild(species);
% annotation.setAttribute('annotation');
