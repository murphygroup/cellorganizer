function G_psf = setsingMTinten_hela(G)

mtt = 19;
psfpeak = 0.08;

G_psf = G*(mtt/psfpeak);

end
