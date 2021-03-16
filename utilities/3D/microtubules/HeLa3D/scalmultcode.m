function scalmult = scalmultcode()
%What does this code do? It doesn't seem to be called anywhere - Devin Sullivan: 2/15/12

cc = zeros(11,11,11);
cc(:,6,6) = 1;
cc2 = psf_blur_hela_mean(cc);
psfpeak = max(cc2(:));

scalmult = avgmtint/psfpeak;
