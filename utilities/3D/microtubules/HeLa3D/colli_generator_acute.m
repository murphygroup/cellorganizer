function [G,tcheck] = colli_generator_acute(n,mu_len,sigma_len,colli_min_number,Dboth,imgcent_coordinate,XYZres,centroradius,stepsize,rand_seed,cellnum,subfolder,batchno)

% Dec 28, 2010: This code quits if the next point within the collinearity cone cannot be generated even after 100 tries.

% tcheck returns true if can be generated, return false if cannot be generated

tcheck = false; imXYZ = []; mtXYZ = []; G = []; randlengths = [];  %%
save(['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '-tcheck.mat'], 'tcheck', 'imgcent_coordinate', 'imXYZ', 'mtXYZ', 'G', 'n','randlengths');  %%


tcheck = true;
colli_array = [-1:0.05:0.05,collinearityrange(20,0.1)];
colli_array = [-1:0.05:1];

% reset random number generator
rand('seed',rand_seed);
randn('seed',rand_seed);

% microtubule coordinates
G = zeros(size(Dboth,1),size(Dboth,2),size(Dboth,3));

if subfolder==1  %%
% generate the lengths based on truncated normal distribution
phi1 = normcdf_my((0-mu_len)/sigma_len);
phi2 = 1;

randlengths = mu_len + sigma_len*(sqrt(2)*erfinv(2*(phi1+(phi2-phi1)*rand(n,1))-1));
end

if subfolder==2  %%
% generate the lengths based on Erlang distribution
randlengths = randraw('erlang',[stepsize ceil(mu_len/(stepsize))],[n 1]);  %% subfolder = 2;
end


randlengths = ceil(randlengths./XYZres(1));

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
[maxstorelength,maxstoreidx] = max(lens);
store_idx = 1;

% generate the first two points for all the microtubules
for j = 1:n
  if randlengths(j) ~= 0
    mcheck = 0;
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
	while 1
        	thet = acos(((2*rand) - 1));
        	phi = 2*pi*rand;
        	% sample points from a pdf that looks quadratic
        	centroradii = ((rand)^(1/3))*centroradius;
        	newpoint = [centroradii*sin(thet)*cos(phi);centroradii*sin(thet)*sin(phi);centroradii*cos(thet)]; 
        	tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres);
		if  ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) ... 
				&& (tempmt2(2) >= 1)  && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))
			if (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)
				mtXYZ{j}(:,1) = newpoint;
          			imXYZ{j}(:,1) = tempmt2;
				t = 1;
				if randlengths(j) == 1
					break;
				end
					
          			while 1 % This generates the second point
%            thet = pi*rand;
            				thet = acos(((2*rand) - 1));
	    				phi = 2*pi*rand;
            				rand_vect = [stepsize*sin(thet)*cos(phi);stepsize*sin(thet)*sin(phi);stepsize*cos(thet)];
             				newpoint = mtXYZ{j}(:,1) + rand_vect;
                  
             				% convert the random vector to image coordinate
             				tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres);
					if  ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) ...
							&& (tempmt2(2) >= 1)  && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))

		   	         		if (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)
							mtXYZ{j}(:,2) = newpoint;
                					imXYZ{j}(:,2) = tempmt2;
							t = 2;
							break;
             					end % End of If
					end % End of If
          			end % End of local while
          			break; % break out from outermost while.
       			end % End of If
		end % End of If
	end % End of while for first two points


	colli_min = colli_min_number;
	minnum = stepsize*colli_min;

	% high_bound sets the bound to the boundary of the nuclear, plasma membrane
	% the value is in microns

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
            		tempmt2 = cart2img(newpoint, imgcent_coordinate, XYZres);
      
            		% check if this coordinate is within bounds
	    		% this already satisfies the collinearity parameter since the point was generated to satisfy the collinearity parameter
	    		if ((tempmt2(1) <= size(Dboth,1)) && (tempmt2(2) <= size(Dboth,2)) && (tempmt2(1) >= 1) && (tempmt2(2) >= 1) && (tempmt2(3) >= 1) && (tempmt2(3) <= size(Dboth,3)))
	    			% Sep 10, 2009: Above if statement added because I have no padding along X, Y
            			if (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 0)
                			mtXYZ{j}(:,t) = newpoint;
                			imXYZ{j}(:,t) = tempmt2;
					colli_min = colli_min_number;
                                        minnum = stepsize*colli_min;
					break;
            			elseif (Dboth(tempmt2(1),tempmt2(2),tempmt2(3)) == 1) && colli_min < 0.55
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
			 if m > 100
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
        if mcheck > 50
		tcheck = false;
		G = [];
        	return;
	end
    end % End of while
  end % end of if
end % end of biggest for .. for j = 1:n


% Now we generate G
b = sum(abs(convn(imXYZ{1},[1 -1])),1);
b(end) = [];
fin_lens = [];
for j = 1 : n
	b = sum(abs(convn(imXYZ{j},[1 -1])),1);
	b(end) = [];
        t = 1;
        G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) = G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) + 1;
        fin_lens(j) = t;
        for t = 2:size(imXYZ{j},2)
		if b(t) ~= 0
			G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) = G(imXYZ{j}(1,t),imXYZ{j}(2,t),imXYZ{j}(3,t)) + 1;
                        fin_lens(j) = t;
                end
        end
end


save(['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '.mat'], 'imgcent_coordinate', 'imXYZ', 'mtXYZ', 'G', 'n','randlengths');

delete(['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mu_len) '_siglen_' num2str(sigma_len) '_colli_' num2str(colli_min_number) '-tcheck.mat']);  %%

end
