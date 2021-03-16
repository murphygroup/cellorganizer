function [obj,protimg] = tp_gengaussobjs( varargin )
% Generate a set of gaussian objects from the model
% Graphical Model: Y <-- X --> Z
%This method can be run in two modes. The mode is determined by the number
%of input arguments.
%
%If the variables protimg and newpos ARE NOT specified the function will
%generate gaussian discs saved as a cell array within the output 'obj'.
%These discs are sampled from the model mean and the generated cov. below.
%They are cut off at 2 standard deviations. protimgsampled will be returned
%[]
%
%If the variables protimg and newpos ARE specified, the function will
%return an image containing the sampled gaussians in the 'protimgsampled'
%image. 'obj' and img will be returned []

% Author: Tao Peng
%
% Copyright (C) 2008-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 9, 2012 I. Cao-Berg Changed report of synthesized objects
% March ??, 2012 D. Sullivan
%   - added sampling method
%   - added direct object addition
%   - added trimmed method
% March ??, 2012 I. Cao-Berg Added varargin structure
% March 13, 2012 D. Sullivan
%   - added user specified sampling density as varargin{6}
%   - fixed trimmed contrast stretching
%   - fixed trimmed point mass
%   - fixed trimmed object 'add' instead of 'replace'
% March 25, 2012 D. Sullivan
%   - changed sampling method to be uniform random within vesicle rather
%   than gaussian
%   - restructured code to optimize efficiency
% August 4, 2012 D. Sullivan Removed "trimmed" method, no longer supported, fixed
%                            memory bugs
% August 4, 2012 D. Sullivan Moved object sampling inside for loop so that
%                            the objects don't have to be held in
%                            memory before being added
% August 4, 2012 D. Sullivan Added check that a sampled object is
%                            representative of the distribution
% August 7, 2012 D. Sullivan Reverted back to adding objects to the images
%                            after they were all sampled.
%
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

i = 0;

if(nargin<2)
    %icaoberg 8/6/2012
    error('Wrong number of input arguments');
end

model = varargin{1};
n = varargin{2};

protimg = varargin{3};
newpos = varargin{4};
objectmethod = varargin{5};
%D. Sullivan 3/14/13 added param structure to list of varargin. need to
%call it before setting the param.samplingDensity
param = varargin{7};
param.samplingDensity = varargin{6};

try
    debug = param.debug;
catch
    debug = false;
end

%icaoberg 8/6/2012
clear varargin;

if ~(strcmpi(objectmethod,'sampled')|| ...
        strcmpi(objectmethod,'disc'))
    warning('CellOgranizer: Unrecognized object fill method, defaulting to disc');
    objectmethod = 'disc';
end

obj = [];
img = [];

%fprintf( 1, '%s%%', '0');

%devins 8/6/2012
%D. Sullivan 3/14/13 changed to normal param structure
% paramobj.imagesize = size(protimg);
param.imagesize = size(protimg);
invalid_objs = 0;

% devins 8/7/2012
% while i < n
%changed from i<n to i<=n
while i <= n
    % Generate covariance matrices
    PofSizex=0;
    PofSizey=0;
    PofSizez=0;
    objtry = -1;
    while ( PofSizex < 0.1 || PofSizex > 0.9 || PofSizey < 0.1 ||...
            PofSizey > 0.9 ||PofSizez < 0.1 || PofSizez > 0.9 )
        objtry = objtry+1;
        
        x = ml_rnd(model.x);
        %x = 2;
        model.y_x.mu = model.y_x.a1*(1-exp(-model.y_x.b1*x));
        model.z_x.mu = model.z_x.a1*(1-exp(-model.z_x.b1*x));
        model.y_x.sigma = model.y_x.a2*(1-exp(-model.y_x.b2*x));
        model.z_x.sigma = model.z_x.a2*(1-exp(-model.z_x.b2*x));
        y = ml_rnd(model.y_x);
        z = ml_rnd(model.z_x);
        
        %         if y>x || z>x || z>y || x<=0 || y<=0 || z<=0
        %             continue
        %         end
        
        %     newSigma = diag([x y z].^2)/25;
        newSigma = diag([x y z].^2);
        
        %if any of the sigmas are less than 1 pixel, set them to 1.
        %this prevents point masses
        if(sum(newSigma(newSigma>0)<1)>0)
            newSigma(newSigma~=0&newSigma<1) = 1;
        end
        
        PofSizex = logncdf(x,model.x.mu,model.x.sigma);
        PofSizey = normcdf(y,model.y_x.mu,model.y_x.sigma);
        PofSizez = normcdf(z,model.z_x.mu,model.z_x.sigma);
        
       

% grj 6/30/13 comment out unused code. 
%         mumat = newpos(i,:);
%         
%         %need to check that the object is the appropriate size
%         %D. Sullivan 3/14/13 seems silly to generate at 6sigma if we only need 2.
%         %Further, if more sigma are desired, we will get the wrong answer
%         % imageSize = ceil([6*sqrt(sigma(1,1)),...
%         %     6*sqrt(sigma(2,2)),6*sqrt(sigma(3,3))]);
%         imageSize = ceil([6*sqrt(newCov(1,1)),...
%             6*sqrt(newCov(2,2)),6*sqrt(newCov(3,3))]);

        
    end
    
    %grj 6/30/13 move rotation code out of build-a-valid-gaussian loop.
    % Random rotation
    t = rand * pi;
    Rx = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
    t = rand * pi;
    Ry = [cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];
    t = rand * pi;
    Rz = [cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
    R = Rz*Ry*Rx;
    newCov = R*newSigma*R';
    if i ==0
        i = 1;
    end
    
    invalid_objs = invalid_objs+objtry;
    %generate gaussian object (vesicle)
%     obj{i} = ml_gaussobj(newCov,param);
    
    obj{i} = ml_gaussimg(newCov, param)
    
    %      obj = ml_gaussobj(newCov,paramobj);
    if (strcmpi(objectmethod,'sampled'))
        numsamples{i} = ceil(param.samplingDensity*sum(abs(newCov(:))));%c*inten_coeff;
        %randomly sample the multivariate gaussian
    end
    
    i = i+1;
    
    % devins 8/7/2012 reverted back to adding objects after they have all been sampled
    % for speed now that the memory issues have been solved.
    %      if(strcmpi(objectmethod,'disc'))
    %             protimg = ml_imaddobj2(protimg,obj,...
    %             struct('method','replace','pos',newpos,'objectmethod',objectmethod));
    %      elseif(strcmpi(objectmethod,'sampled'))
    %
    %             %randomly sample the multivariate gaussian
    %             numsamples = ceil(param.samplingDensity*sum(abs(newCov(:))));%c*inten_coeff;
    %
    %             %properly fill the vesicles based on the method
    %             paramsampling = struct('method','add','pos',newpos(i,:),...
    %             'objectmethod',objectmethod);
    %             paramsampling.numsamplescell = numsamples;
    %             protimg = ml_imaddobj2(protimg,obj,paramsampling);
    %
    %      end
    
    %indicate an object has been successfully added
    %fprintf( 1, [ repmat('\b', 1, length(num2str(100*(i-1)/n))+1) '%s'], ...
    %    repmat(' ', 1, length(num2str(100*(i-1)/n))+1) );
    %fprintf( 1, [ repmat('\b', 1, length(num2str(100*(i-1)/n))+1) '%s%%'], num2str(100*i/n) );

    %fprintf( 1, [ repmat('\b', 1, length(num2str(floor(100*(i-1)/n))+1)) '%s'], ...
    %    repmat(' ', 1, length(num2str( floor(100*(i-1)/n))+1) ) );

    %fprintf( 1, [ repmat('\b', 1, length(num2str( floor(100*(i-1)/n))+1) ) '%s%%'], num2str(floor(100*i/n)) );
    %fprintf( 1, '%s', '.');
end
%fprintf( 1, '%s\n', '' );

%devins 8/7/2012 reverted back to adding objects after they were sampled
if(strcmpi(objectmethod,'disc'))
    protimg = ml_imaddobj2(protimg,obj,...
        struct('method','replace','pos',newpos,'objectmethod',objectmethod));
elseif(strcmpi(objectmethod,'sampled'))    
    %properly fill the vesicles based on the method
    paramsampling = struct('method','add','pos',newpos,...
        'objectmethod',objectmethod);
    paramsampling.numsamplescell = numsamples;
    protimg = ml_imaddobj2(protimg,obj,paramsampling);
    
end

%devins 8/6/2012
if debug
    disp(['The number of invalid objects was: ' num2str(invalid_objs)]);
end
%

%do the contrast stretching before converting to uint8
protimgtmp = ml_bcimg(protimg,[],[0 255]);
protimg = uint8(protimgtmp);
