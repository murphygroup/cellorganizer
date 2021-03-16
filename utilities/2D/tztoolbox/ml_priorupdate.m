function [posterior, class, iteration] = ml_priorupdate(p, O, alpha_values, distancematrix)
% [posterior, class, iteration] = ML_PRIORUPDATE(p, O, alpha_values, distancematrix)
% This function implements the prior updating method, which is the heart of
% this work.  Given the evidence of each cell from the single-cell SVM
% classifier, for each cell the algorithm updates the likelihoods by
% considering the prior distributions of each class and generates the
% posterior probabilty.
%
% Input:
%   p is the evidence matrix of dimension
%	{ number_of_classes X number_of_cells }
%   O is the label vector obtained from the single-cell classifier
%   alpha_values is the strength parameter which controls the extent to
%   which the prior probabilities are adjusted at each iteration 
%   distancematrix is a matrix of dimension
%	{ number_of_cells X number_of_cells }
%	showing the connectivity of the cells (which nodes in the
%	graph are connected); 1 implies connected
%
% Output:
%   posterior is the posterior probability distribution of the class
%	for each cell
%   class is the updated labels 
%   iteration is the number of iterations
%
%   Usage example:     ml_priorupdate(p, O, alpha_values, distancematrix)

% version 1, February 25, 2006
%
% This file is part of the software package described in "A graphical model
% approach to automated classification of protein subcellular location
% patterns in multi-cell images" by Shann-Ching Chen and Robert F. Murphy,
% BMC Bioinformatics 7:90 (2006).
%
% Copyright (C) 2006  Shann-Ching Chen and Robert F. Murphy,
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

[N, ClassNum] = size(p);    
total_test_num = N;

flag = ones(total_test_num ,1);
conf = [];
for J = 1:total_test_num
    p_prior_incorporated = p;
    for I = 1:total_test_num                  
        if flag(I) == 0
            p_prior_incorporated(I,:) = 0;
        end
    end
    
    for K = 1:total_test_num
        second_max = zeros(total_test_num ,1);                     
        for I = 1:total_test_num
            [Y,i] = max(p_prior_incorporated(I,:));
            p_temp = p_prior_incorporated(I,:);
            p_temp(i) = -inf;
            [Y,i] = max(p_temp);
            second_max(I) = Y;
        end                     
    end
    
    p_prior_incorporated = p_prior_incorporated - second_max*ones(1,ClassNum);
    [Y,i] = max(p_prior_incorporated');
    [most_confident,index_c] = max(Y);
    
    if most_confident == 0
        mess_cell = find(flag == 1);
        index_c = mess_cell(1);
        i = O(index_c)*ones(1,total_test_num);          
        
        % others are set to zero
        setzero = [1:ClassNum];
        setzero(O(index_c)) = [];                  
        p(index_c,setzero) = 0;%%
    end
    
    class = i(index_c);
    flag(index_c) = 0;
    
    conf = [conf ;[most_confident,index_c,class]];
end %  for J = 1:total_test_num

flag = ones(total_test_num ,1);
class_previous = O';            
class = zeros(total_test_num ,1);%                        

conf = sortrows(conf,2);
iteration = 0;   
while 1                     
    iteration = iteration + 1;                                   
    con = conf(:,2);
    for J1 = 1:total_test_num
        if isempty(find(distancematrix(:,con(J1))==1))
            % classify it without any prior information
            [Y,class(con(J1))] = max(p(con(J1),:));
            flag(class(con(J1))) = 0;
            if Y == 0
                class(con(J1)) = O(con(J1));
            end   
        else
            cc = find(distancematrix(:,con(J1))==1);
            position = [];
            for J2 = 1:total_test_num
                position = [position ; find(cc==con(J2))]; 
            end      
            c1 = cc(position);
            
            allconf = 0;
            for J2 = 1:size(c1,1)
                allconf = allconf + conf(find(conf(:,2) == c1(J2)),1); 
            end            
            
            priorconf = zeros(1,ClassNum);
            for J2 = 1:size(c1,1)
                if class_previous(c1(J2)) < ClassNum + 0.5 
                    priorconf(class_previous(c1(J2))) = priorconf(class_previous(c1(J2))) + conf(find(conf(:,2) == c1(J2)),1);            
                end
            end
            
            prior = (1/ClassNum)*ones(1,ClassNum); 
            for J2 = 1:ClassNum
                if priorconf(J2) > 0
                    prior_new(J2) = (prior(J2) + priorconf(J2)*alpha_values)/(1 + allconf*alpha_values);      
                else
                    prior_new(J2) = prior(J2)/(1+allconf*alpha_values);             
                end
            end   
            prior = prior_new; 
            p_prior_incorporated = (ones(total_test_num,1)*prior).*p;      
            [Y,class(con(J1))] = max(p_prior_incorporated(con(J1),:));  
            flag(class(con(J1))) = 0;  
            if Y == 0
                class(con(J1)) = class_previous(con(J1));
            end
        end                              
    end
    
    if sum(class_previous ~= class) == 0
        %					[begin_prior;prior]
        break;               
    else
        f1 = find(class_previous-class ~=0);
        for f3 = 1:length(f1)
            f2 = find(conf(:,2) == f1(f3));
            conf(f2,1) = 0;	% if the label is changed, the confidence should be set to zero      
            setzero = [1:ClassNum];
            setzero(class(f1(f3))) = [];                  
            p(f1(f3),setzero) = 0;%%                  
        end
        class_previous = class;
    end                            
end %while 1

posterior = p_prior_incorporated;
