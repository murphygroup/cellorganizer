function [idx_sda, ignoreidx] = ml_stepdisc(features, logfilename, Debug)
% [idx_sda, ignoreidx] = ml_stepdisc(features)
%
% Perform Stepwise Discriminant Analysis (SDA) on FEATURES
% where features is a cell array of feature matrices,
% one cell per class. It returns a set of ranked informative features
% in idx_sda. The SDA log file is stored in 'sda.log'
% Reference:  Jennrich, R.I. (1977), "Stepwise Discriminant Analysis,"
% in Statistical Methods for Digital Computers, eds. 
% K. Enslein, A. Ralston, and H. Wilf, New York: John Wiley & Sons, Inc.
% Input:
% 	features: a cell array of length ClassNum
%		in which each cell contains a feature matrix of
%		dimensions imagenums X number of features
% 	ClassNum: the number of classes
%	classnames: a cell array of length ClassNum containing the
%		names of the classes 
% 	imagenums: a vector of length ClassNum containing the
%		number of images in each class
% Output:
%   idx_sda: the ranked feature index
%   ignoreidx: the feature index which are constant or weird.  These
%       features are ignored for SDA analysis

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

%   30-Aug-2006 Written by S. C. Chen

idx_sda = [];
if ~iscell(features)
    fprintf(1,'Please arrange features into a cell array\n');
    return;    
end

% When Debug is 0, only the sda index and the sda.log will be generated
% When Debug is 1, in addition to calculate the index and the log file, 'Statistics
% for Removal' and 'Statistics for Entry' will also be displayed on the
% standard output.
if nargin == 1
    logfilename = 'sda.log';
    Debug = 0;
end

if nargin == 2    
    Debug = 0;
end

if nargin == 3
    if Debug ~= 1
        Debug = 0;
    end
end

q = size(features,2);       % Number of class (groups)
m = size(features{1},2);    % Total number of features (variables)
ng = zeros(1,q);            % Number of features in each class
for i=1:q
    ng(i) = size(features{i},1);
    features{i} = double(features{i});  % change the features in double format to increase the numerical precision
end

% The following Boxes and Circles are part of the flowchart in the reference.

% Box 1~8: initialization, only have to read the data once
% Box 1
xgibar = zeros(q,m);
wij = zeros(m,m);
lambda = 1;
% Box 3,4,5,6,7
for g=1:q 
    for t=1:ng(g)
        % Box 2
        xgti = features{g}(t,:);
        Di = xgti - xgibar(g,:);
        xgibar(g,:) = ((t-1) * xgibar(g,:) + xgti)/t;
        Ei = xgti - xgibar(g,:);
        for i=1:m
            for j=i:m
                wij(i,j) = wij(i,j) + Di(i)*Ei(j);
            end % j=i:m
        end % i=1:m
    end % t=1:ng(g)
end % g=1:q

% Box 8
n = sum(ng);
xibar = zeros(1,m);
for g=1:q
    xibar = xibar + xgibar(g,:)*ng(g);
end
xibar = xibar/n;

tij = wij;
for i=1:m
    for j=i:m
        for g=1:q            
            tij(i,j) = tij(i,j) + ng(g) * (xgibar(g,i) - xibar(i)) * (xgibar(g,j) - xibar(j));
        end % g=1:q
    end % j=i:m
end % i=1:m
% Go to Circle 1

% Circle 1
% Box 9
for i=1:m
    for j=i:m
        wij(j,i) = wij(i,j);        tij(j,i) = tij(i,j);        
    end % j=i:m    
end % i=1:m
Di = diag(wij)';
Ei = diag(tij)';

ignoreidx=find(Di==0);
% if any idx in the Di is 0, its feature is wierd  and
% should be ignored in the feature selection
if ~isempty(ignoreidx)
    feat = features;
    for i=1:size(features,2)
        feat{i}(:,ignoreidx) = [];
    end   
    idx_sda = ml_stepdisc(feat, logfilename, Debug);   
    for i=1:length(ignoreidx)
        idx_sda(idx_sda>ignoreidx(i)) = idx_sda(idx_sda>ignoreidx(i)) + 1;
    end
    return;
end

% Tolerance threshold, default 0 in SAS
TOL = 0;
% Significance Level to Enter, default 0.15 in SAS
SigIn = 0.15;
% Significance Level to Stay, default 0.15 in SAS
SigOut = 0.15;

iteration = 1;  % counter for interation

% Log file for iteration number, FeatNumber, FValue, ProbabilityF, and WilksLambda
flog = fopen(logfilename,'w');

p = 0;  
while 1
    Circle = 1;    
    % Box 10
    FLAG = 0;
    divtmp = diag(wij)./(Di');    
    if p == 0
        % Box 14
        if isempty(find(divtmp >= TOL))   
            % Go to Circle 2
            Circle = 2;
        else
            % Box 15     
            divtmp2 = diag(wij)./diag(tij);                             
            [Vk, k] = min(divtmp2(find(divtmp>=TOL)));
            k = find(divtmp2 == Vk);
            
            df1 = q-1;  df2 = n-p-q+1;    F_dist_value = (1-Vk)*(n-p-q)/(Vk*(q-1));   PrF = 1-fcdf(F_dist_value,df1,df2);     
            % Box 16           
            if PrF <= SigIn
                % Box 17         
                FLAG = -1;
                p = p + 1;
                % Go to Circle 3
                Circle = 3;
            else
                % Go to Circle 2
                Circle = 2;
            end
        end
    else
        % Box 11
        divtmp2 = diag(wij)./diag(tij);                   
%        divtmp2(find(diag(wij) >= 0)) = MAXINT;
        [Vk, k] = min(divtmp2(find(diag(wij)<0)));

        k = find(divtmp2 == Vk);          
        
        df1 = q-1;  df2 = n-p-q+1;    F_dist_value = (Vk-1)*(n-p-q+1)/(q-1);   PrF = 1-fcdf(F_dist_value,df1,df2);                 
        if PrF > SigOut
            % Box 13         
            FLAG = 1;
            p = p - 1;
            % Go to Circle 3
            Circle = 3;       
        else
            % Box 14            
            if isempty(find(divtmp >= TOL))   
                % Go to Circle 2
                Circle = 2;
            else
                % Box 15     
                divtmp2 = diag(wij)./diag(tij);                                 
                [Vk, k] = min(divtmp2(find(divtmp>=TOL)));
                k = find(divtmp2 == Vk);
            
%                 divtmp2(find(divtmp < TOL)) = MAXINT;
%                 [Vk, k] = min(divtmp2);
                
                % Box 16      
                df1 = q-1;  df2 = n-p-q;    F_dist_value = (1-Vk)*(n-p-q)/(Vk*(q-1));   PrF = 1-fcdf(F_dist_value,df1,df2);    
                if PrF <= SigIn              
                    % Box 17  
                    FLAG = -1;
                    p = p + 1;
                    % Go to Circle 3
                    Circle = 3;
                else
                    % Go to Circle 2
                    Circle = 2;
                end  
            end       
        end
    end % end if p == 0
    if Circle == 2
        % Box 21~26,  doing classification & stop iteration
        break;
    end
    if Circle == 3
        % Box 18
        %tz+ 28-Oct-2006
        k = k(1);
        %tz++
        
        wij = SWP(wij, k, FLAG);
        tij = SWP(tij, k, FLAG);
        lambda = Vk * lambda;
        
        % Box 19
        F_to_remove = zeros(1,m);
        F_to_enter = zeros(1,m);
        
        fprintf(flog,'%3d\t',iteration);     iteration = iteration + 1;      
        
        if Debug == 1        
            fprintf(1,'iteration %3d:  ',iteration-1);         
            fprintf(1,'\n\nStatistics for Entry, DF = %d, %d\n\n', df1, df2);
        end            
        for i=1:m
            if wij(i,i) > 0
                F_to_enter(i) = (tij(i,i) - wij(i,i)) * (n-p-q)/(wij(i,i)*(q-1));     
                if Debug == 1
                    fprintf(1,'VAR%3d\t\t%.2f\n',i+1,F_to_enter(i));    %%%                     
                end
            end
        end
        if Debug == 1        
            fprintf(1,'\n\n');
        end  
        
        if Debug == 1        
            fprintf(1,'Statistics for Removal, DF = %d, %d\n\n', df1, df2);
        end        
        for i=1:m
            if wij(i,i) <= 0
                F_to_remove(i) = (wij(i,i) - tij(i,i)) * (n-p-q+1)/(tij(i,i)*(q-1));
                if Debug == 1                
                    fprintf(1,'VAR%3d\t%.2f\n',i+1,F_to_remove(i));    %%%    
                end
            end
        end               
        if Debug == 1        
            fprintf(1,'\n\n');
        end  
        
        % PrF is to test the significance of the change in lambda resulting from the addition of u                
        if PrF < 0.0001
            PrFstr = '<.0001';
        else
            PrFstr = sprintf('%.4f',PrF);
        end
        
        if Debug == 1              
            if FLAG == -1
                % entere a varaible;   k+1 is compatible for SAS's format
                fprintf(1,'VAR%3d (F=%.2f) is entered\n',k+1,F_to_remove(k));    
            else
                % remove a varaible;   k+1 is compatible for SAS's format                          
                fprintf(1,'==>VAR%3d (F=%.2f) is removed\n',k+1,F_dist_value);           
            end        
        end
        
        if FLAG == -1
            % entere a varaible;   k+1 is compatible for SAS's format
            fprintf(flog,'%3d\t\t%.2f\t%s\t%.8f\n',k+1,F_to_remove(k),PrFstr,lambda);    
        else
            % remove a varaible;   k+1 is compatible for SAS's format          
            fprintf(flog,'\t%3d\t%.2f\t%s\t%.8f\n',k+1,F_dist_value,PrFstr,lambda);           
        end        
    end
end % end while 1

fclose(flog);

% Read the data into a matrix
[ Step, FeatNumber, FValue, ProbabilityF, WilksLambda ] = textread(logfilename , '%d %d %f %s %f');
% adjust the k+1 index to k index
FeatNumber = FeatNumber - 1;           
% Select the "Good" features as being the ones having Less than 
% 0.0001 probability of getting an F value of FValue or larger
% if the null hypthesis (that the classes are the same) was true.
% I.e. we just go down the list to the place where the first one
% with ProbabilityF greater than <.0001 is.
good_feature_idx = [];
for i=1:length(ProbabilityF)
    PrF = ProbabilityF{i};
    switch(PrF)
        case '<.0001'
	        good_feature_idx = [good_feature_idx i];
        otherwise
    end
end
best_subset = FeatNumber(good_feature_idx);
% Finally, sort these "best" features by their F-value or Wilks' lambda
fval=FValue(good_feature_idx);
[SortedF,idx]=sort(fval);
idx_sda=best_subset(flipdim(idx,1))';


function swpmat = SWP(mat, k, flag)
% The sweep operation, a special type of exchange operator similar to those
% used in linear programming, matrix inversion, and elsewhere.
% Reference:  Jennrich, R.I. (1977), "Stepwise Regression,"
% in Statistical Methods for Digital Computers, eds. 
% K. Enslein, A. Ralston, and H. Wilf, New York: John Wiley & Sons, Inc.

[row,col] = size(mat);
swpmat = zeros(row,col);

if flag == -1
    % forward sweep
    swpmat(k,k) = -1/mat(k,k);
    for i=1:row
        if i~=k
            swpmat(i,k) = mat(i,k)/mat(k,k);
        end
    end
    for j=1:col
        if j~=k
            swpmat(k,j) = mat(k,j)/mat(k,k);
        end
    end    
    for i=1:row
        for j=1:col
            if i==k || j==k
                continue;
            else
                swpmat(i,j) = mat(i,j)-mat(i,k)*mat(k,j)/mat(k,k);
            end
        end
    end
else 
    if flag == 1
        % inverse sweep
        swpmat(k,k) = -1/mat(k,k);
        for i=1:row
            if i~=k
                swpmat(i,k) = -mat(i,k)/mat(k,k);
            end
        end
        for j=1:col
            if j~=k
                swpmat(k,j) = -mat(k,j)/mat(k,k);
            end
        end    
        for i=1:row
            for j=1:col
                if i==k || j==k
                    continue;
                else
                    swpmat(i,j) = mat(i,j)-mat(i,k)*mat(k,j)/mat(k,k);
                end
            end
        end        
    end
end
