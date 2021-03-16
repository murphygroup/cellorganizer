function redcdpars = reducedlibrary(cellnum,subfolder)

%subfolder = 1;

load intensities mcXs
I = 1;
for cellno = [1 2 4 6 10 12:42 46:49 51 52]
	mcX = mcXs(cellno,:); 
	tt(I)= mcX(4); 
	I = I + 1;
end
mtt = mean(tt);
stt = std(tt);

% Compute total intensity per point
q = zeros(11,11,11);
q(5,5,5) = 1;
q = psf_blur_hela_mean(q);
q = setsingMTinten_hela(q);

% From real image limits
ll = mtt - stt;
ul = mtt + stt;

protim3 = getrealimage_hela(cellnum);
UL = sum(protim3(:))/((sum(q(:))/mtt)*ll); % Total tubulin in units of that used while generating images
LL = sum(protim3(:))/((sum(q(:))/mtt)*ul);

batchno = 1; % for getting the synthetic image randlengths

%n_inputArray = [5,50:50:350];
%mulen_inputArray  = [5 10 15 20 25 30 35];
n_inputArray = [5, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350];
mulen_inputArray  = [5 10 15 20 25 30 35 40];
colli_inputArray = [0.9 0.95 0.98];

if subfolder==1
coeff_var_Array = [0 0.1 0.2 0.3];
else
coeff_var_Array = [0];
end

I = 1;
for n = n_inputArray
	for mulen = mulen_inputArray
		for coeffvar = coeff_var_Array
			sigmalen = coeffvar*mulen;
			for colli = colli_inputArray
				myfile =  ['outputs_' num2str(subfolder) '/images/cell_' num2str(cellnum) '/batch_' num2str(batchno) '/new_n_' num2str(n) '/new_sim_n_' num2str(n) '_mulen_' num2str(mulen) '_siglen_' num2str(sigmalen) '_colli_' num2str(colli) '.mat'];
				if exist(myfile,'file') == 2
					[randlengths] = getrandlengths_hela(n,mulen,sigmalen,colli,cellnum,batchno,subfolder);
					sp = sum(randlengths);
                                
					if (sp > LL) && (sp < UL)
						redcdpars(I,:) = [n mulen sigmalen colli];
						I = I + 1;
					end
				end
			end
		end
	end
end

end % Function ends
