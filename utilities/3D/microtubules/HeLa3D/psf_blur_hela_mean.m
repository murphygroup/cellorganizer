function G_psf = psf_blur_hela_mean(G)

load HPA_285

G_psf = imfilter(G,PSF,'conv','same');

end
