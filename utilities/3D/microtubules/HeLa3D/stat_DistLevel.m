function feat_vector = stat_DistLevel(image,imnum)

[levels] = getDistLevel2(image,imnum);

% Below is modified from MI_histprops.m
%--------------------------------------------------------------------------
% 1. Calculate histogram 
%--------------------------------------------------------------------------
NL = length(levels);
% Generate gray level vector
% intialize parameters
Gray_vector = zeros(1,NL);
Histogram = zeros(1,NL);
for i =1:NL
    %Gray_vector(i) = sum(levels{i}(:,4));  %%~Gray_vector~Normalized
    Gray_vector(i) = nanmean(levels{i}(:,4),1);
    %Histogram(i) = size(levels{i},1);
    Histogram(i) = 1;  %%
end
%Gray_vector = Gray_vector/(sum(Gray_vector));  %%~Gray_vector~Normalized

%--------------------------------------------------------------------------
% 2. Now calculate its histogram statistics
%--------------------------------------------------------------------------
% Calculate obtains the approximate probability density of occurrence of the intensity
% levels
Prob                = Histogram./(nansum(Histogram));
% 2.1 Mean 
Mean                = nansum(Prob.*Gray_vector);
% 2.2 Variance
Variance            = nansum(Prob.*(Gray_vector-Mean).^2);
% 2.3 Skewness
Skewness            = calculateSkewness(Gray_vector,Prob,Mean,Variance);
% 2.4 Kurtosis
Kurtosis            = calculateKurtosis(Gray_vector,Prob,Mean,Variance);

%-------------------------------------------------------------------------
% 3. Insert all features and return
%--------------------------------------------------------------------------
feat_vector =[Mean Variance Skewness Kurtosis];
% End of funtion


%--------------------------------------------------------------------------
% Utility functions
%--------------------------------------------------------------------------
function Skewness = calculateSkewness(Gray_vector,Prob,Mean,Variance)
% Calculate Skewness
term1    = Prob.*(Gray_vector-Mean).^3;
term2    = sqrt(Variance);
Skewness = term2^(-3)*sum(term1);

function Kurtosis = calculateKurtosis(Gray_vector,Prob,Mean,Variance)
% Calculate Kurtosis
term1    = Prob.*(Gray_vector-Mean).^4;
term2    = sqrt(Variance);
Kurtosis = term2^(-4)*sum(term1);

function [I, nl, gl] = ParseInputs(varargin)
% parsing parameter checking
% Inputs must be max seven item
iptchecknargin(1,5,nargin,mfilename);
%
% Check I
I = varargin{1};
iptcheckinput(I,{'logical','numeric'},{'2d','real','nonsparse'}, ...
    mfilename,'I',1);
% ------------------------
% Assign Defaults
% -------------------------
%
if islogical(I)
    nl = 2;
else
    nl = 256; % Modified by Aabid
end
gl = getrangefromclass(I);

% Parse Input Arguments
if nargin ~= 1

    paramStrings = {'NumLevels','GrayLimits'};

    for k = 2:2:nargin

        param = lower(varargin{k});
        inputStr = iptcheckstrs(param, paramStrings, mfilename, 'PARAM', k);
        idx = k + 1;  %Advance index to the VALUE portion of the input.
        if idx > nargin
            eid = sprintf('Images:%s:missingParameterValue', mfilename);
            msg = sprintf('Parameter ''%s'' must be followed by a value.', inputStr);
            error(eid,'%s', msg);
        end

        switch (inputStr)
            case 'NumLevels'
                nl = varargin{idx};
                iptcheckinput(nl,{'logical','numeric'},...
                    {'real','integer','nonnegative','nonempty','nonsparse'},...
                    mfilename, 'NL', idx);
                if numel(nl) > 1
                    eid = sprintf('Images:%s:invalidNumLevels',mfilename);
                    msg = 'NL cannot contain more than one element.';
                    error(eid,'%s',msg);
                elseif islogical(I) && nl ~= 2
                    eid = sprintf('Images:%s:invalidNumLevelsForBinary',mfilename);
                    msg = 'NL must be two for a binary image.';
                    error(eid,'%s',msg);
                end
                nl = double(nl);

            case 'GrayLimits'

                gl = varargin{idx};
                iptcheckinput(gl,{'logical','numeric'},{'vector','real'},...
                    mfilename, 'GL', idx);
                if isempty(gl)
                    gl = [min(I(:)) max(I(:))];
                elseif numel(gl) ~= 2
                    eid = sprintf('Images:%s:invalidGrayLimitsSize',mfilename);
                    msg = 'GL must be a two-element vector.';
                    error(eid,'%s',msg);
                end
                gl = double(gl);
        end
    end
end

