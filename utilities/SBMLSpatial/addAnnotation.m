function annotation = addAnnotation(docNode)

annotation = docNode.createElement('annotation');
cellorg = docNode.createElement('cellorganizer');
annotation.appendChild(cellorg);
% annotation.setAttribute('annotation');
