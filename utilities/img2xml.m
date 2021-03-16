function img2xml( img , filename )
%INDEX2XML Converts an indexed image to a sparse matrix, and saves as an xml
%file
%
%Input arguments    Description
%img                image, represented as a matrix, to save as an xml
%filename           Location to save .xml file to

% Author: Michelle Mackie (mmackie@andrew.cmu.edu)
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

if( nargin > 2)
    error('CellOrganizer: Wrong number of input arguments.' );
end

if (length(size(img)) ~= 3)
   error('Invalid image' ); 
end

docNode = com.mathworks.xml.XMLUtils.createDocument('sparse');


for z_index=1:size(img,3)
    z_slice = docNode.createElement('z-slice');
    img_slice = img(:,:,z_index);
    img_sparse = sparse(img_slice);
    [i,j,val] = find(img_sparse);

    for n=1:length(val)
       x = i(n);
       y = j(n);
       temp_val = val(n);

       sparse_element = docNode.createElement('element');
       x_node = docNode.createElement('x');
       x_text = docNode.createTextNode(num2str(x));
       x_node.appendChild(x_text);
       sparse_element.appendChild(x_node);

      y_node = docNode.createElement('y');
      y_text = docNode.createTextNode(num2str(y));
      y_node.appendChild(y_text);
      sparse_element.appendChild(y_node);

      val_node = docNode.createElement('value');
      val_text = docNode.createTextNode(num2str(temp_val));
      val_node.appendChild(val_text);
      sparse_element.appendChild(val_node);
      
      z_slice.appendChild(sparse_element);
      
    end
    docNode.getDocumentElement.appendChild(z_slice);
end

xmlwrite([filename '.xml'], docNode);
end
