function post=ml_label2post(label,N)

%ML_LABEL2POST converts class integer labels to nx1 coding
%   POST=ML_LABEL2POST(LABEL) converts the vector LABEL 
%   to a matrix POST. The ith row and kth column of POST is
%   1 and all other elements in its ith row is 0 if the ith 
%   element of LABEL equals to k. If a row has all 0s, it will
%   be converted into 0.
%
%   POST=ML_LABEL2POST(LABEL,N) specifies the the number categories by N.
%   So POST is a nxN matrix if LABEL is a nx1 vector. 
%   
%   Example: If X =[1  ml_label2post(x,3) is [1 0 0  
%                   2                         0 1 0            
%                   0]                        0 0 0]             
%                                           
%   See also ML_POST2LABEL

% Copyright (C) 2006  Murphy Lab
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

%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments
%   10-FEB-2005 Modified TINGZ
%       - change algorithm to increase the speed
%   16-Mar-2006 Modified T. Zhao
%       - support 0 label

%tz- 16-Mar-2006 
% if any(label==0)
%     error('label could not be zero')
% end
%tz--

if length(size(label))>3 | (size(label,1)>1 & size(label,2)>1)
    error('label should be a vector')
end

% tic
% post=zeros(length(label),max(label));
% 
% for i=1:length(label)
%     post(i,label(i))=1;
% end
% toc
maxlabel=max(label);
id=eye(maxlabel);

if nargin==2  
    %tz+ 16-Mar-2006 
    if N>size(id,2)
    %tz++
        id=[id,zeros(size(id,1),N-size(id,2))];
    %tz+ 16-Mar-2006    
    else
        if N<size(id,2)
            warning('The labels exceed specified number of categroies');
        end
    end
    %tz++
end   

%Add a row of zeros at the beginning for 0 labels
id = [zeros(1,size(id,2));id];

post=id(label+1,:);