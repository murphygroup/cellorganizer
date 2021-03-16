function [marginals] = ml_priorupdatevp(ev, edges, lambda)
% [marginals] = ML_PRIORUPDATEVP(ev, edges, lambda)
% This function implements the prior updating method with the voting potential.
% Given the evidence of each cell from the single-cell SVM
% classifier, for each cell the algorithm updates the likelihoods by
% considering the prior distributions of each class and generates the
% posterior marginal probabilty.
%
% Input:
%   ev is the evidence matrix of dimension
%	{ number_of_classes X number_of_cells }
%   edges is a matrix of dimension
%	{ number_of_cells X number_of_cells }
%	showing the connectivity of the cells (which nodes in the
%	graph are connected); 1 means connected
%   lambda is the smoothing parameter which controls the how much we want
%   to belief from the prior distribution and the information from the
%   neighbor nodes; the smaller the lambda, the larger we want to believe
%   the neighbor information
%
% Output:
%   marginals is the posterior probability distribution of the class
%	for each cell
%
%   Usage example:     ml_priorupdatevp(ev, edges, lambda)
% 25 Mar 2006
% Written by Shann-Ching Chen

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

% Copyright (C) 2006  Shann-Ching Chen and Robert F. Murphy,
% Carnegie Mellon University

[N, nclass] = size(ev);    
message = ones(nclass, N, N)/nclass;
phy = ones(nclass, N)/nclass;
belief = zeros(nclass, N);  

while 1
    messagenew = message;
    for X = 1:N             
        % Mxx: message from clique to center node x
        neighbor = find(edges(:,X)==1);
        m = length(neighbor);
        if m == 0
            continue;
        else              
            message_sum = zeros(nclass,1);
            for Vi = 1:m 
                message_sum = message_sum + message(:,neighbor(Vi),X);
            end
            phy(:,X) = (message_sum + lambda/nclass)/(m+lambda);
            phy(:,X) = phy(:,X)/sum(phy(:,X));

            phi = ev(X,:)' .* phy(:,X); 
            phi = phi/sum(phi);
            
            % update the messages from node i to other cliques
            for Vi = 1:m
                message(:,X, neighbor(Vi)) = phi;
            end            
            
        end
    end
    
    m = abs(message - messagenew);
    if sum(m(:)) < 0.00001
        break;
    else
        messagenew = message;
    end        
end    

for X = 1:N
    neighbor = find(edges(:,X)==1);
    m = length(neighbor);
    if m == 0
        belief(:,X) = ev(X,:)';
    else              
        belief(:,X) = ev(X,:)' .* phy(:,X);        
    end
    belief(:,X) = belief(:,X)/sum(belief(:,X));
end

marginals = belief';
