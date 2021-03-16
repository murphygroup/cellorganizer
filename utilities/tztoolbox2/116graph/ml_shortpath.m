function nodes = ml_shortpath(gm,param)
%ML_SHORTPATH Find shortest path in a directed acyclic graph
%   NODES = ML_SHORTPATH(GM) returns the shortest path from the first node to
%   the last node in the binary graphical matrix GM. The searching is 
%   accelerated by dynamic programming.
%   
%   NODES = ML_SHORTPATH(GM,PARAM) specifies the parameters by the structure
%   PARAM, which has the following fields:
%       'gtype' - type of graphical matrix
%           1 : binary (default)
%           2 : continuous weights
%   
%   See also

%   30-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('gtype',1));
nnode = size(gm,2);
nodes = [];

switch param.gtype
    case 1 %binary distance
        gm(end) = 1;

        for i=nnode:-1:1
            %check columns
            pos = find(gm(:,i)>0);
            if length(pos)>1
                gm(pos,i) = gm(i,i);
            end

            %check rows
            for j=1:length(pos)
                d = gm(pos(j),i)+1;
                if gm(pos(j),pos(j))==0 | gm(pos(j),pos(j))>d
                    gm(pos(j),pos(j)) = d;
                end
            end
        end

        nodes = 1;
        pos = 1;
        while nodes(end)<nnode
            pos = find(gm(pos,:)==gm(nodes(end),pos)-1);
            nodes = [nodes pos(1)];
%             pos = max(find(gm(:,pos(1))==gm(nodes(end-1),nodes(end))));
            pos = pos(1);
        end
    case 2 %continuous distance
        for i=nnode:-1:1
            pos = find(gm(1:i-1,i)>0);
            if length(pos)>=1
                gm(pos,i) = gm(i,i)+gm(pos,i);
            end
            %check rows
            for j=1:length(pos)
                d = gm(pos(j),i);
                if gm(pos(j),pos(j))==0 | gm(pos(j),pos(j))>d
                    gm(pos(j),pos(j)) = d;
                end
            end
        end
        nodes = 1;
        pos = 1;
        while nodes(end)<nnode
            pos = find(gm(pos,pos+1:end)==gm(pos,pos))+pos;
            nodes = [nodes pos(1)];
            pos = pos(1);
        end

end

