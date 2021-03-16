function [G,imXYZ,mtXYZ,randlengths,resolution] = ... 
     MT_synth(numberMT,mu_len,colli_min_number,...
     Dboth,imgcent_coordinate,XYZres,method, rand_seed, param)
%Inputs
%------
%n = number of microtubules
%mu_len = mean length of microtubules
%colli_min_number = collinearity of microtubule segments
%Dboth = cellular outline (negative mask e.g. cell area = 0 and exterior = 1
%   Dboth corresponds to the addition of the cell_image logical and the
%   dna_image logical however is itself a double (due to the addition)
%imgcent_coordinate = I'm guessing this is the center of the image/cell?
%XYZres = assuming this is the image resolution in microns (0.2)?
%method = the method you wish to use. If this is left blank the full bounce
%method will be used. 
%method = 1 --> Aabid Shariff Gaussian 
%method = 2 --> Erlang, one step rebound 
%method = 3 --> Erlang, multi-step rebound 
%
%Outputs
%-------
%G = the raw microtubule image
%tcheck = removed by DPS 2/17/12
%imXYZ =
%mtXYZ =
%randlengths 

% Author: Jieyue Li
%   - originally called: "colli_generator_acute_jl2"
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% February 15, 2012 Devin Sullivan
%   - changed input and output arguments 
%   - made usable for synthetic images 
%   - commented parameters (see below)
% March 9, 2012 Devin Sullivan
%   - added downsampling
%   - added upsampling
%   - added removal of added rows/columns
% March 12, 2012 I. Cao-Berg Updated license and documentation
% March 18, 2012 Devin Sullivan
%   - added contrast stretching
% April 12, 2012 Devin Sullivan
%   - added old (Aabid Shariff) method
%   - added method option 
% July 25, 2012 R.F. Murphy
%   - increase contrast stretching (centrosome too bright)
% July 26, 2012 D. Sullivan
%   - updated method to perform bilinear interpolation at image resize
% Feb. 5-6 2013 Devin Sullivan 
%   - modified code to work for the multiple resolutions in the future
%   rather than hardcoding the up and down sampling. Additional
%   modifications will be needed to allow variables to work at different
%   resolutions or pass model specific variables
% Feb 27, 2013 D. Sullivan 
%   Changed the way resolutions are dealt with. Now dealt with in
%   model2instance.m
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

% Dec 28, 2010: This code quits if the next point within the collinearity cone cannot be generated even after 100 tries.
%adaptation of colli_generator_acute_jl2 by Jieyue Li
%Header added 2/15/12 Devin Sullivan


%previous input parameters that are now fixed or were not explicitly passed
%to the function by Jieyue

%If no random seed is specified, get a random random-seed
if ~exist('rand_seed', 'var') | isempty(rand_seed)
    rng('default')
    rng('shuffle')
else
    rand('seed',rand_seed);
    randn('seed',rand_seed);
end

if ~exist('param', 'var') | isempty(param)
    param = struct();
end


default_param = struct( ...
    'quit_if_cant_synth', true, ...
    'maxiter', 100);
%eventually move in centroradius, stepsize

param = ml_initparam(param, default_param);

%fixed by DPS 2/17/12
%centroradius = 0.4 hard coded by Jieyue
centroradius = 0.4;
stepsize = 0.2;% hard coded by Jieyue
%rand_seed = a number to fix the random seed with DPS 2/15/12 fixed at 5
% rand_seed = 5;
%cellnum = a number used to indicate what cell we are on for saving. removed by DPS 2/15/12
%subfolder = a number indicating which output version we are saving to. removed by DPS 2/15/12
%batchno = 1 hard coded reference by Jieyue. removed by DPS 2/15/12

%D. Sullivan 2/5/13 add check for resolution, set to default if not
%included
if ~isfield(XYZres,'objects')
    warning('no object resolution specified for microtubule model. Defaulting to [0.2,0.2,0.2]')
    XYZres.objects = [0.2,0.2,0.2];
end
if(~exist('method','var'))
    method = 'bounce_control';
end

%D. Sullivan 2/27/13 keep track of object synthesis
%Note: removed all resolution adjustment code
resolution = XYZres.objects;


%turn it back into a logical 
Dboth = AdjustResolutions(Dboth,XYZres.cell,XYZres.objects,1,'round');
imgcent_coordinate = floor(imgcent_coordinate.*XYZres.centrosome./XYZres.objects);
% Dboth = logical(Dboth);
%End 3/9/12 addition
%%

%not explicitly passed to the function by Jieyue. the if ~exist statements
%removed by DPS 2/17/12
%sigma_len = std of microtubules (Seems to only be used for file naming since Jieyue's update,
%removed by DPS 2/18/12)
%   nosave = 0;
%    nomcheck = 0; %removed by DPS 2/17/12
colli_thre = 0.34;
%  angle_thre = 100;  %% degree
angle_thre = 111;  %% degree

%D.Sullvian 2/5/13 make sure the object model is trained on the same
%resolution in every dimension for the checkback to work 
if exist('XYZres.objects')&&~all(XYZres.objects==XYZres.objects(1))
    error('Inconsistent resolutions used for MT model. Currently CellOrganizer only supports equal resolutions for XYZ');
end
%   checkback = round(6/XYZres(1)); % = 30; stepsize = XYZres(1) = 0.2, 6 (microns)/0.2 = 30 (pixels or steps).
%checkback = ceil(min(size(Dboth,1),size(Dboth,2))/256*6/XYZres.objects(1)); %
checkback = 6/XYZres.objects(1); %

% tcheck returns true if can be generated, return false if cannot be generated

%tcheck = false;
imXYZ = []; mtXYZ = []; G = []; randlengths = [];  %%


%tcheck = true;
%colli_array = [-1:0.05:0.05,collinearityrange(20,0.1)];%what is this line for, immediately rewritten, commenting out DPS 2/15/12
colli_array = [-1:0.05:1];

% reset random number generator
%rand('seed',rand_seed);
%randn('seed',rand_seed);

% microtubule coordinates
%Dboth gives coordinates of the cell shape:(comment inserted DPS 2/15/12)
G = zeros(size(Dboth,1),size(Dboth,2),size(Dboth,3));


% if n = 0, no MT in cell, so all other values = 0(comment inserted DPS 2/15/12)
if numberMT == 0
    mu_len = 0;
    mtXYZ = G;
    imXYZ = mtXYZ;
    randlengths = 0;
    return
end

%%%%
%temporary test DPS 3/9/12
%mu_len = 2*mu_len;
%%%%

%figure out what method we are using and initialize the random lengths 
switch method
    case 'no_bounce_control'
        % generate the lengths based on Erlang distribution
        randlengths = randraw('erlang',[stepsize ceil(mu_len/(stepsize))],[floor(numberMT) 1]);
    case 'bounce_control'
        % generate the lengths based on Erlang distribution
        %randlengths = randraw('erlang',[stepsize ceil(mu_len/(stepsize))],[n 1]);  %% subfolder = 2; [nx1] array
        randlengths = randraw('erlang',[stepsize ceil(mu_len/(stepsize))],[floor(numberMT) 1]);

    otherwise
        warning(['No known method bounce control selected, defaulting to multistep bounce control']);
        randlengths = randraw('erlang',[stepsize ceil(mu_len/(stepsize))],[floor(numberMT) 1]);
end
%no longer using the Gaussian model
% if method==1  %%
% % generate the lengths based on truncated normal distribution
% %so since the current model does not take a variance parameter for length,
% %we will set the coeffvar arbitrarily to 0.1 which is within the range
% %[0,0.1,0.2,0.3] swept by the original code 
% %DPS 4/12/12
% coeffvar = 0.3;
% sigma_len = coeffvar*mu_len;
% phi1 = normcdf_my((0-mu_len)/sigma_len);
% phi2 = 1;
% 
% randlengths = mu_len + sigma_len*(sqrt(2)*erfinv(2*(phi1+(phi2-phi1)*rand(numberMT,1))-1));


randlengths = ceil(randlengths./XYZres.objects(1));

randlengths = sort(randlengths,1,'descend');

% if (min(randlengths) <= 2)
%	error('Length of each MT must be atleast three points')
% end

% Initialize mtXYZ and imXYZ sizes
for j = 1:size(randlengths,1)
    mtXYZ{j} = zeros(3,randlengths(j));
    imXYZ{j} = zeros(3,randlengths(j));
end

ims =[];
mts = [];
lens = 0;
[maxstorelength,maxstoreidx] = max(lens);%this is zero given previous line: DPS 2/15/12
store_idx = 1;

% generate the first two points for all the microtubules
disp(['Generating first two points for all ',num2str(numberMT),' microtubules']);
for j = 1:numberMT
    fprintf('.');
    if mod(j,50)==0
        fprintf('\n');
    end
    if randlengths(j) ~= 0
        mcheck = 0;
        counter1=0;
        while 1 % this loop will end when we are happy about the length generated.
            if randlengths(j) <= maxstorelength
                imXYZ{j} = ims{maxstoreidx}(:,1:randlengths(j));
                mtXYZ{j} = mts{maxstoreidx}(:,1:randlengths(j));
                ims{maxstoreidx} = [];
                mts{maxstoreidx} = [];
                lens(maxstoreidx) = 0;
                [maxstorelength,maxstoreidx] = max(lens);
                break;
            end

            % generate the first two points for all the microtubules
            counter2=0;
            while 1
                %random angles at which to grow? DPS 2/15/12
                thet = acos(((2*rand) - 1));
                phi = 2*pi*rand;
                % sample points from a pdf that looks quadratic
                centroradii = ((rand)^(1/3))*centroradius;
                newpoint = [centroradii*sin(thet)*cos(phi);centroradii*sin(thet)*sin(phi);centroradii*cos(thet)];
                tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres.objects(1));

                if  ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) ...
                        && (tempmt2(2) >= 1)  && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))
                    if (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)
                        mtXYZ{j}(:,1) = newpoint;
                        imXYZ{j}(:,1) = tempmt2;
                        t = 1;

                        if randlengths(j) == 1
                            break;
                        end
                        
                        counter3=0;
                        while 1 % This generates the second point
                            %			              thet = pi*rand;
                            thet = acos(((2*rand) - 1));
                            phi = 2*pi*rand;
                            rand_vect = [stepsize*sin(thet)*cos(phi);stepsize*sin(thet)*sin(phi);stepsize*cos(thet)];
                            newpoint = mtXYZ{j}(:,1) + rand_vect;
                            
                            % convert the random vector to image coordinate
                            tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres.objects(1));
                            if  ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) ...
                                    && (tempmt2(2) >= 1)  && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))
                                
                                if (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)%if there is not already a MT in the way: (comment added DPS 2/15/12)
                                    mtXYZ{j}(:,2) = newpoint;
                                    imXYZ{j}(:,2) = tempmt2;
                                    t = 2;
                                    break;
                                end % End of If
                            end % End of If
                            counter3=counter3+1;
                        end % End of local while
                        break; % break out from outermost while.
                    end % End of If
                end % End of If
                counter2 = counter2+1;
            end % End of while for first two points            
            
            colli_min = colli_min_number;
            minnum = stepsize*colli_min;
            
            % high_bound sets the bound to the boundary of the nuclear, plasma membrane
            % the value is in microns
            
            %Now we generate the remainder of the microtubule
            breakoutflag = false;
            if randlengths(j) > t
                for t = 3:randlengths(j)
                    m=0;
                    KLKL = 0;
                    while 1 % loop until we get the right coordinates
                        m = m+1;
                        % generate a random vector
                        %            		thet = pi*rand;
                        thet = real(acos((((stepsize - minnum)*rand) + minnum)/stepsize));
                        phi = 2*pi*rand;
                        rand_vect = [stepsize*sin(thet)*cos(phi);stepsize*sin(thet)*sin(phi);stepsize*cos(thet)];
                        % Transform the random vector such that it lies along the growing microtubule.
                        
                        T = mtXYZ{j}(:,t-1) - mtXYZ{j}(:,t-2);
                        Y = [0;0;stepsize];
                        
                        % we have to rotate Y to T
                        % we use Rodrigues' rotation formula
                        
                        rotangle = acos((Y'*T)/(stepsize^2));
                        
                        temp = cross(Y,T);
                        rotaxis = temp./norm(temp);
                        
                        Amat = [0, -rotaxis(3), rotaxis(2);
                            rotaxis(3),0,-rotaxis(1);
                            -rotaxis(2),rotaxis(1),0];
                        
                        rotmat = eye(3) + (Amat.*sin(rotangle)) + ((Amat^2).*(1-cos(rotangle)));
                        
                        newpoint = mtXYZ{j}(:,t-1) + (rotmat*([rand_vect(1),rand_vect(2),rand_vect(3)]'));
                        
                        % convert the random vector to image coordinate
                        tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres.objects(1));
                        
                        % check if this coordinate is within bounds
                        % this already satisfies the collinearity parameter since the point was generated to satisfy the collinearity parameter
                        if ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) && (tempmt2(2) >= 1) && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))
                            % Sep 10, 2009: Above if statement added because I have no padding along X, Y
                            if (strcmpi(method,'no_bounce_control') && Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)
                                mtXYZ{j}(:,t) = newpoint;
                                imXYZ{j}(:,t) = tempmt2;
                                colli_min = colli_min_number;
                                minnum = stepsize*colli_min;
                                break;
                            elseif (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0) && checkRebounce(tempmt2,imXYZ{j}(:,max(1,t-1-checkback+1):(t-1)),angle_thre)  %%
                                mtXYZ{j}(:,t) = newpoint;
                                imXYZ{j}(:,t) = tempmt2;
                                colli_min = colli_min_number;
                                minnum = stepsize*colli_min;
                                break;
                                %elseif (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 1) && colli_min < colli_thre
                            elseif colli_min < colli_thre  %%
                                % if it hits edge, store the info, and do not add it to the MTs.
                                lens(store_idx) = t-1;
                                mts{store_idx} = mtXYZ{j}(:,1:t-1);
                                ims{store_idx} = imXYZ{j}(:,1:t-1);
                                % mtXYZ{j} = []; % i commented this, because the new set will be overwritten
                                % imXYZ{j} = [];
                                [maxstorelength,maxstoreidx] = max(lens);
                                store_idx = store_idx + 1;
                                breakoutflag = true;
                                mcheck = mcheck + 1;
                                m = 0;
                                break; % breaks from the while loop to generate until right point
                            else
                                if tempmt2(3) ~= imXYZ{j}(3,t-1)
                                    KLKL = KLKL + 1;
                                    colli_min = colli_array(end-KLKL+1);
                                    minnum = stepsize*colli_min;
                                    m=0;
                                end
                            end
                        end
                        if m > 100  %%
                            KLKL = KLKL + 1;
                            colli_min = colli_array(end-KLKL+1);
                            minnum = stepsize*colli_min;
                            m=0;
                        end
                        
                    end % End of while
                    if breakoutflag
                        break;
                    end
                end % End of For .. for t = 3: ..
            end % end of if
            if (t == randlengths(j)) && (~breakoutflag)
                break;
            end
            
            %Commented this back in and increased iterations. GRJ 11/9/14
            if mcheck > param.maxiter        %%
                tcheck = false;  %%
                if param.quit_if_cant_synth % Causes program to stop, why????! DPS 2/18/12, Re above: The quit is because there may be no guaranteed solution. GRJ 11/9/14
                    G = [];
                    return;
                else
                    break
                end
            end
            counter1 = counter1+1;
        end % End of while
    end % end of if
end % end of biggest for .. for j = 1:n


% Now we generate G
b = sum(abs(convn(imXYZ{1},[1 -1])),1);
b(end) = [];
fin_lens = [];
fprintf('\nGenerating microtubules\n');
for j = 1 : numberMT
    fprintf('.');
    if mod(j,50)==0
        fprintf('\n');
    end
    b = sum(abs(convn(imXYZ{j},[1 -1])),1);
    b(end) = [];
    t = 1;
    G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) = G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) + 1;
    fin_lens(j) = t;
    for t = 2:size(imXYZ{j},2)
        if (b(t) ~= 0) && (prod([imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)])~=0)
            G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) = G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) + 1;
            fin_lens(j) = t;
        end
    end
end
fprintf('.\n');

%D. Sullivan 2/27/13 Removed resizing code (now done in model2instance.m)


%%
%Added DPS 3/18/12
%contrast stretching
low_high = [0 max(G(:))];
G = ml_bcimg(G,low_high,[0 255]);
%end 3/18/12 contrast stretching addition
%%
end
