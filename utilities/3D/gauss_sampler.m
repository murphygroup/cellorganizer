function img = gauss_sampler(mu,sigma,img,param)
% GAUSS_SAMPLER Samples from a given gaussian object to give sampled values to the gaussian 

% Author: Devin Sullivan (devins@cmu.edu)
% March 12, 2012 I. Cao-Berg Added license and modified documentation
% March 12, 2012 I. Cao-Berg Added param.sampling.density
% March 14, 2012 Devin Sullivan Fixed sampling density
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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


%DPS 2/29/12

%determine how many samples to take. 
%this is proportional to the inten_coeff
%c = ??;%factor by which to multiply intensity coefficient 
% 
% try
%    samplingDensity = param.samplingDensity;
% catch
%    samplingDensity = 100;
% end
for k = 1:length(param.numsamplescell)
  try
      numsamples = param.numsamplescell{k};
  catch
      numsamples = ceil(samplingDensity*sum(abs(sigma(:))));%c*inten_coeff;
  end
  %numsamples = ceil(100*sum(abs(sigma(:))));%c*inten_coeff;

  %randomly sample the multivariate gaussian
  %mu is mean, sigma is covariance 
  %this should return a multivariate gaussian approximating 
  try 
      Rnumrnd =param.Rnumrndcell{k};
  catch
       Rnums = mvnrnd(mu,sigma,numsamples);
      Rnumrnd = ceil(Rnums);
  end


  %make sure we have 3D
  if(size(Rnumrnd,2)~=3)
      return
  end

  count = 0;
  %given the postitions of the samples, add a 1 to each sample point within
  %the gaussian object
  for i = 1:numsamples 
      flag = 0;
      %need to make sure that there are no zero values in Rnumrnd
      for j = 1:size(Rnumrnd,2)

          if(Rnumrnd(i,j)<=0||Rnumrnd(i,j)>size(img,j))
              %disp('out of bounds value suggested, skipping sample');
              flag = 1;
          end
      end
      if(flag == 1)
          continue
      end
        
      img(Rnumrnd(i,1),Rnumrnd(i,2),Rnumrnd(i,3)) = img(Rnumrnd(i,1),Rnumrnd(i,2),Rnumrnd(i,3))+1;
  end

end
