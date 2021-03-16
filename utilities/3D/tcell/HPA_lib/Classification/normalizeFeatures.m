function data = normalizeFeatures( data);
% NORMDATA = FEATNORM( DATA)
% 
% Normalizes data by dividing row elements by the magnitude of the 
% respecitve rows
% 
% Last modified Justin Newberg 11 February 2009

data = data./repmat( sqrt(sum(data(:,:).^2,2)), [1 size(data,2)]);
