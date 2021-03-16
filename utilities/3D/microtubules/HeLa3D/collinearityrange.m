function colli_array = collinearityrange(numofpoints,finalcolli) 


% clear all
%close all

% numofpoints = 5;
% stepsize = 0.09;

% finalcolli = 0.97;

% NOTE: I am using an approximate cone. It does not have a curved surface.

% This is volume of cone in terms of height and theta

% finalvolume = (1/3)*pi*((stepsize)^3)*((tan(acos(finalcolli)))^2);

% We want volumes of cone (for 5 pts) to be [0, (2^0)*V, (2^1)*V, (2^2)*V, (2^3)*V], where V = initial volume
% Initial volume = finalvolume/2^(numofpoints-3);

% But we do not need to compute any volume, since we only want collinearity. We can set the condition that volume reduces by half every time.
colli_array = zeros(1,numofpoints);

colli_array(1) = finalcolli;

for I = 2:(numofpoints-1)
	colli_array(I) = cos(atan((sqrt(0.5))*(tan(acos(colli_array(I-1))))));
end
% initialvolume = finalvolume/(2^(numofpoints-2));

colli_array(numofpoints) = 1;
