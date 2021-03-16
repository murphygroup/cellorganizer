%cooc3d
%Authors: Carl Philips and Daniel Li (2008)
%http://facweb.cs.depaul.edu/research/vc/contact.htm
%
%Modified by Gregory R. Johnson 2/18/14
%   gj@andrew.cmu.edu
%   Mutltiple image inputs for co-occurance across images
%
%Syntax
%[featureVector, coocMat] = cooc3d (data, 'distance', [1;2;4;8], ...
%   'direction', [0 1 0;1 1 0; 0 1 -1])
%
%Description
%reads in a vector of cubes and outputs a matrix of haralick features and 
%the 3D Co-Occurrence matrices.
%
%Input:
%Data: a vector of cubes with the fourth dimension identifying the cube.
%data(:,:,:,1) = rand(20,20,20); %cube 1
%data(:,:,:,2) = rand(20,20,20); %cube 2
%
%Parameters:
%numGray: Integer indicating the number of graylevels to use when
%performing the graylevel resizing (rescaling).
%distance: a nx1 array of distances that will be used when analyzing the
%image. Default is [1,2,4,8];
%direction: a nx3 array of direction offsets in [row, column, vertical]
%format. The vertical value increases from top to bottom 
%   the standard 2D directions 
%   [0 1 0]    0 degrees
%   [-1 1 0]   45 degrees
%   [-1 0 0]   90 degrees
%   [-1 -1 0]  135 degrees
%
%   The additional 9 directions that make this a 3D Co-Occurrence ...
%       algorithm
%             horizontal, vertical
%   [0 1 -1]   0 degrees, 45 degrees
%   [0 0 -1]   straight up
%   [0 -1 -1]  0 degrees, 135 degrees
%   [-1 0 -1]  90 degrees, 45 degrees
%   [1 0 -1]   90 degrees, 135 degrees
%   [-1 1 -1]  45 degrees, 45 degrees
%   [1 -1 -1]  45 degrees, 135 degrees
%   [-1 -1 -1] 135 degrees, 45 degrees
%   [1 1 -1]   135 degrees, 135 degrees
%   Default is all 13 directions.
%
%Output:
%featureVector = haralick values for each cube (this is what's used for
%                classification. Each row pertains to a different cube.
%coocMat = the Co-Occurrence matrices.
%          coocMat(y,x,direction,distance,cube number)

%featureVector(:,1:12) = the haralick features for distance 1, ...
%   direction 1;
%featureVector(:,13:24) = the haralick features for distance 1, ...
%   direction 2;
%featureVector(:,127:167) = the haralick features for distance 2, ...
%   direction 1;
%The haralick features used (in order) are:
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, Inertia, 
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability, 
%Inverse Variance.


%Designed and tested for cubes with axis 20 voxels long.


%function [featureVector,coocMat] = cooc3d (data,distance, directions)
function [featureVector,coocMat] = cooc3d (im1, im2, distance, numLevels)
%inputStr = {'Distance','Direction'};

%Default settings
if ~exist('distance', 'var') | isempty(distance)
    distance = [1,2,4,8,16]; %more or fewer distances?
end

if ~exist('numLevels', 'var') | isempty(numLevels)
    numLevels = 64;
end

numHarFeature = 12; %changing this may break the harFeatures function


offsets{1} = [0  1  0;  %xy plane 
             -1  1  0; 
             -1  0  0; 
             -1 -1  0];          

if ndims(im1) > 2

    offsets{2} = [0  0 -1]; %z axis

    % offsets{3} = [0 -1 -1; ... ... %diagonal up directions
    %              -1  0 -1; ...
    %               1  0 -1; ...
    %              -1  1 -1; ...
    %               1 -1 -1; ...
    %              -1 -1 -1; 
    %               1  1 -1];

    offsets{3} = [1  0 -1; ...
                  1 -1 -1; ...
                  0 -1 -1; ...
                 -1 -1 -1; ...
                  1  0  1; ...
                  1 -1  1; ...
                  0 -1  1; ...
                 -1 -1  1];
end

if ~exist('im2', 'var') | isempty(im2)
    im2 = im1;
end


numoffsets = length(offsets); %number of directions, currently 13
coocMat = zeros(numLevels, numLevels, numoffsets, size(distance,2));

% for iteration=1:size(data,4) %each new cube

tempVec = zeros(1,0);
for i = 1:length(offsets)
    
    for dist =1:length(distance) %distance
        
        [harMat, coocMat(:,:,i,dist)] = graycooc3d(...
            im1, im2, distance(dist),numLevels,...
            numHarFeature,offsets{i}); 
        temphar = zeros(1,0);
        
        %organizing the data so each cube's data is on a row
        for clicks =1:size(harMat,1)
            temphar = cat(2,temphar,harMat(clicks,:));
        end
        tempVec = cat(2,tempVec,temphar);
              
    end
end
    %produces a larger space to separate cubes
    %haralickMat = cat(1,haralickMat,space3);
    featureVector = tempVec;

%     
%     disp(['completed cube number' num2str(iteration)])
% end
return


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%I = the 3D image matrix
%distance = a vector of the distances to analyze in
%numLevels = the number of graylevels to be used
%numHarFeature = the number of haralick features to compute
%offSet = a matrix of the directions to analyze in
%
%harMat = a matrix of the haralick features in the format harMat(direction,
%feature)
%coMat the Co-Occurrence matrices produced
function [harMat,coMat]= graycooc3d(im1, im2,distance,numLevels,numHarFeature,...
    offSet)

%**************Variable initialization/Declaration**********************
harMat =0;


noDirections = size(offSet,1); %number of directions, currently 13
coMat = zeros(numLevels,numLevels,noDirections);


%************************graylevel resizing*******************************
numLevels = numLevels-1; %don't touch. Logical adding issue.
minImage = min(vertcat(im1(:), im2(:)));
im1=im1-(minImage);
im2=im2-(minImage);
% min(min(min(im1)));
maxImage = max(vertcat(im1(:), im2(:)));
tempShift = double(maxImage)/double(numLevels);
im1 = floor(double(im1)/double(tempShift));
im2 = floor(double(im2)/double(tempShift));
im1=im1+1;
im2=im2+1;

numLevels = numLevels+1; %don't touch. Logical adding issue.
if max(max(max(im1))) > numLevels
    disp('Error is graylevel resizing.')
    disp('cooc3d.m');
    return
end


%**************************Beginning analysis*************************
%Order of loops: Direction, slice, graylevel, graylevel locations
for direction =1:noDirections %currently 13 (for the 3d image)

    tempMat = zeros(numLevels,numLevels,size(im1,3));
    for slicej =1:size(im1,3)
         for j=1:numLevels %graylevel
             
             %finds all the instances of that graylevel
            [rowj,colj] = find(im1(:,:,slicej)==j);  

            %populating the Cooc matrix.
            for tempCount = 1:size(rowj,1) 
                rowT = rowj(tempCount) + distance*offSet(direction,1);
                colT = colj(tempCount) + distance*offSet(direction,2);
                sliceT = slicej+distance*offSet(direction,3);
                [I1, I2, I3] = size(im1);
                if rowT <= I1 && colT <= I2 && sliceT <= I3
                    if rowT > 0 && colT > 0 && sliceT > 0
                        
                        %Error checking for NANs and Infinite numbers
                        IIntensity = im2(rowT,colT,sliceT);
                        if ~isnan(IIntensity)
                            if ~isinf(IIntensity)
                                %Matlab doesn't have a ++ operator.
                                tempMat(j,IIntensity,slicej)= tempMat...
                                    (j,IIntensity,slicej)+1;
                            end
                        end
                    end
                end
            end
         end
    end
    for slicej =1:size(im1,3)
        coMat(:,:,direction)= coMat(:,:,direction)+tempMat(:,:,slicej);
    end
end

coMat = sum(coMat,3);
coMat(1) = 0;
%extracting the Haralick features from the Co-Occurrence matrices
harMat = harFeatures(coMat,numHarFeature);
return


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%coMat = Co-occurrence matrices 2D stack upon eachother (i j k) k is number
%of directions analyzed. For 3d that's 13. Created in cooc3d.m
%numHarFeature is the number of variables you will be extracting.
%Haralick order
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, Inertia, 
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability, 
%Inverse Variance.
%harMat = matrix of the haralick features in the format harMat(direction,
%feature)
function [harMat]= harFeatures(coMat, numHarFeature)

%numHarFeature=12;
%numPosFeature=12; %If you add any more features bump this up.
numLevels = size(coMat,1); %number of graylevels
harMat = zeros(numHarFeature,size(coMat,3));
%%%%%%tempHarMat = zeros(numPosFeature,1);  %continue working here....
%tempCoMat=zeros(size(coMat,1),size(coMat,2));


for iteration = 1:size(coMat,3) %directions

    
%%%%%%%%%%%%%%%%%%%%Preparation

%%%%%%%determining various p values

    pij = sum(sum(coMat(:,:,iteration))); %already normalized
    
    if pij == 0
        pij = 1;
    end
    
    coMat(:,:,iteration)=coMat(:,:,iteration)./pij;

    tempmux=0;
    tempmuy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempmux =  tempmux+(i*(coMat(j,i,iteration)));
            tempmuy =  tempmuy+(j*(coMat(j,i,iteration)));
        end
    end
    mux=tempmux; %mux
    muy=tempmuy;

    tempx=0;
    tempy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempx = tempx+ (i-mux)^2*coMat(j,i,iteration);
            tempy = tempy+ (j-muy)^2*coMat(j,i,iteration);
        end
    end
    sigx=tempx; %sigx
    sigy=tempy;
    



%Calculations
    tempEnergy =0;
    tempEntropy=0;
    tempCorr=0;
    tempCont=0;
    tempGen=0;
    tempVar=0;
    tempMean=0;
    tempInert=0;
    tempShade=0;
    tempTen=0;
    tempInVar=0;
    for j=1:numLevels
        for i=1:numLevels
            value = coMat(j,i,iteration);
            
            tempEnergy = tempEnergy+ value^2;
            if(value~=0) 
                tempEntropy = tempEntropy + (value * log10(value));
            end
            tempCorr = tempCorr+ ((i-mux)*(j-muy)*(value/(sigy*sigx)));
            n=(abs(i-j))^2;
            tempCont = tempCont+ value*n;
            tempGen = tempGen+ value/(1+abs(1-j));
            tempVar = tempVar + ((i - mux)^2)*value+((j-muy)^2)*value;
            tempMean = tempMean + (i+j)*(value);
            tempInert = tempInert+ (i-j)^2*(value);
            tempShade=tempShade+ ((i+j-mux-muy)^3)*(value);
            tempTen = tempTen+ (((i + j - mux - muy)^4) .* (value));
            if i~=j
                tempInVar=tempInVar+ value/(i-j)^2;
            end
        end
    end
    harMat(1,iteration)=tempEnergy;         %Energy
    harMat(2,iteration) = -tempEntropy;     %Entropy
    
    if isnan(tempCorr)
        tempCorr = 0;
    end
    
    harMat(3,iteration)=tempCorr;           %Correlation
    harMat(4,iteration)=tempCont;           %Contrast
    harMat(5,iteration) = tempGen;          %Homogeneity
    harMat(6,iteration) = tempVar/2;        %Variance
    harMat(7,iteration)=tempMean/2;         %Sum Mean
    harMat(8,iteration)=tempInert;          %Inertia
    harMat(9,iteration)=tempShade;          %Cluster Shade
    harMat(10,iteration) = tempTen;         %Cluster Tendency
    harMat(11,iteration) = max(max(coMat(:,:,iteration))); %Max Probability
    harMat(12,iteration) = tempInVar;       %Inverse Variance
    
    clear 'tempEnergy' 'tempEntropy' 'tempCorr' 'tempCont' 'tempGen';
    clear 'tempVar' 'tempMean' 'tempInert' 'tempShade';
    clear 'tempTen' 'tempInVar';

end
%makes it so that rows are cases
harMat = harMat';
return
