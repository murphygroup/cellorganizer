function [d, redcov, errortype] = ml_mahal(Y,X)
% ML_MAHAL Mahalanobis distance with correction for zero variance variables.
% 
%   ML_MAHAL(Y,X) gives the Mahalanobis distance of each point in
%   Y from the sample in X (using the MAHAL function from the Matlab
%   STATS toolbox) but with elimination of variables that have zero 
%   variance to avoid inversion problems
%
%   [D, REDCOV, ERRORTYPE] = ML_MAHAL(Y,X) returns the distance in D,
%     the covariance matrix of the reduced variable set in REDCOV, and
%     an error indication in ERRORTYPE (0=no error, NaN= not enough	       %     points to allow Mahal calculation, Inf = at least one variable
%     had zero variance and one point in Y didn't match)
%
%   Created July 20, 2001 by Robert F Murphy
%   Modified July 21, 2001 by Robert F Murphy
%     to handle more cases  
%   Modified January 8, 2002 by Robert F Murphy
%     to trap case where all variables have zero variance
%   Modified January 10, 2002 by Robert F Murphy
%     to fix return of reduced covarance matrix
%   Modified February 26, 2002 by Robert F Murphy
%     to fix case with one variable

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

[rx,cx] = size(X);
[ry,cy] = size(Y);

if cx ~= cy
    error('Requires the inputs to have the same number of columns.');
end

if (rx>1) 
    varX = var(X);
else
    varX = zeros(1,cx);
end;
if (ry>1)
    varY = var(Y);
else
    varY = zeros(1,cy);
end;
redX = [];
redY = [];
j = 0;
for i=1:cx
    % if all of the values for this variable in sample X are the same
    % we can't compute Mahalanobis distance so don't include that variable 
    % in the calculation
    if(varX(i)>0)
        j = j + 1;
        redX(:,j)=X(:,i);
        redY(:,j)=Y(:,i);
    end;
end;
if rx < j
    warning('The number of rows of X must exceed the number of columns with non-zero variance.');
    d = repmat(NaN,ry,1);
    redcov = NaN;
    errortype = NaN;
else
    if (j>0)
        if (j==1)
            redcov=std(redX);
            d=((redY-mean(redX))./redcov(1)).^2;
        else
            d = mahal(redY,redX);
            redcov=cov(redX);
        end
    else
        d = zeros(ry,1);
        redcov = zeros(j);
    end
    % ignoring a zero variance variable will only be accurate if all of the
    % values for that variable in Y are also the same as the value in X,
    % so check to see if any of the points have an infinite distance
    %
    % use a flag so that the warning is only printed once
    errortype = 0;
    for i=1:cx
        if(varX(i)==0 & varY(i)>0)
            if(errortype==0)
                %        warning('Infinite Mahalanobis distance(s) due to zero variance.');
                errortype=Inf;
            end
            for k=1:ry
                % check for points in Y where there is zero variance in a variable in X but
                % where the value for that variable is not the same as the value for X
                if(Y(k,i)~=X(1,i))
                    d(k) = Inf;
                end;
            end;
        end;
    end;
end;



