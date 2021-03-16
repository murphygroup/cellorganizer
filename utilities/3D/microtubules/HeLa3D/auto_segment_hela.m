function auto_segment_hela(cellnum)

cellseg(1,:) = [5 10];
cellseg(2,:) = [6 13];
cellseg(4,:) = [5 14];
cellseg(6,:) = [5 10];
cellseg(10,:) = [6 13];
cellseg(11,:) = [5 11];
cellseg(12,:) = [4 15];
cellseg(13,:) = [6 13];
cellseg(14,:) = [5 9];
cellseg(15,:) = [5 15];
cellseg(16,:) = [7 16];
cellseg(17,:) = [7 13];
cellseg(18,:) = [6 15];
cellseg(19,:) = [5 11];
cellseg(20,:) = [6 12];

cellseg(21,:) = [6 11];
cellseg(22,:) = [6 15];
cellseg(23,:) = [5 14];
cellseg(24,:) = [7 14];
cellseg(25,:) = [5 10];

cellseg(26,:) = [7 15];
cellseg(27,:) = [6 15];
cellseg(28,:) = [5 11];
cellseg(29,:) = [6 15];
cellseg(30,:) = [5 11];
cellseg(31,:) = [5 9];
cellseg(32,:) = [6 11];
cellseg(33,:) = [6 15];
cellseg(34,:) = [6 13];
cellseg(35,:) = [6 15];
cellseg(36,:) = [7 16];
cellseg(37,:) = [6 12];
cellseg(38,:) = [6 12];
cellseg(39,:) = [6 11];
cellseg(40,:) = [4 11];
cellseg(41,:) = [5 11];
cellseg(42,:) = [6 13];
cellseg(46,:) = [6 12];
cellseg(47,:) = [6 11];
cellseg(48,:) = [6 15];
cellseg(49,:) = [7 14];
cellseg(51,:) = [6 11];
cellseg(52,:) = [6 11];


dnaseg(1,:) = [6 9];
dnaseg(2,:) = [7 12];
dnaseg(4,:) = [7 12];
dnaseg(6,:) = [6 9];
dnaseg(10,:) = [7 11];
dnaseg(11,:) = [6 9];
dnaseg(12,:) = [5 11];
dnaseg(13,:) = [7 11];
dnaseg(14,:) = [6 8];
dnaseg(15,:) = [6 14];
dnaseg(16,:) = [8 14];
dnaseg(17,:) = [8 12];
dnaseg(18,:) = [7 14];
dnaseg(19,:) = [6 8];
dnaseg(20,:) = [7 11];

dnaseg(21,:) = [7 10];
dnaseg(22,:) = [7 14];
dnaseg(23,:) = [6 13];
dnaseg(24,:) = [6 13];
dnaseg(25,:) = [6 9];

dnaseg(26,:) = [8 14];
dnaseg(27,:) = [7 14];
dnaseg(28,:) = [6 10];
dnaseg(29,:) = [7 14];
dnaseg(30,:) = [6 11];
dnaseg(31,:) = [6 8];
dnaseg(32,:) = [7 10];
dnaseg(33,:) = [7 13];
dnaseg(34,:) = [7 12];
dnaseg(35,:) = [7 14];
dnaseg(36,:) = [8 15];
dnaseg(37,:) = [7 11];
dnaseg(38,:) = [7 11];
dnaseg(39,:) = [7 10];
dnaseg(40,:) = [5 10];
dnaseg(41,:) = [6 10];
dnaseg(42,:) = [7 12];
dnaseg(46,:) = [7 11];
dnaseg(47,:) = [7 10];
dnaseg(48,:) = [7 14];
dnaseg(49,:) = [8 13];
dnaseg(51,:) = [7 10];
dnaseg(52,:) = [7 10];

% Load images
load(['proc/cell_' num2str(cellnum) '/cell' num2str(cellnum) '.mat'],'cellim3','dnaim3','protim3');

% Load PSFs
infdnaPSF = imfinfo( '3DHeLa_DNA_PSF.tif');
infcellPSF = imfinfo( '3DHeLa_Cell_PSF.tif');

dnaPSF = zeros( infdnaPSF(1).Height, infdnaPSF(1).Width, length(infdnaPSF));
cellPSF = zeros( infcellPSF(1).Height, infcellPSF(1).Width, length(infcellPSF));

for I=1:length(infdnaPSF)
    dnaPSF(:,:,I)=imread('3DHeLa_DNA_PSF.tif',I);
end
for I=1:length(infcellPSF)
    cellPSF(:,:,I)=imread('3DHeLa_Cell_PSF.tif',I);
end

dnaPSF = dnaPSF.^2; % Approximate confocal PSF
cellPSF = cellPSF.^2; 

dnaPSF = ml_downsize(dnaPSF,[4 4 1],'average');
cellPSF = ml_downsize(cellPSF,[4 4 1],'average');
[cell_image,cellPSF2] = deconvblind(cellim3,cellPSF);

[dna_image,dnaPSF2] = deconvblind(dnaim3,dnaPSF);

% Cell image segmentation
segcell = active3Dsegment(cell_image,cellseg(cellnum,1),cellseg(cellnum,2));

% Nucleus image segmentation
segdna = active3Dsegment(dna_image,dnaseg(cellnum,1),dnaseg(cellnum,2));

save(['proc/cell_' num2str(cellnum) '/man_seg.mat'],'segcell','segdna','dnaim3','protim3');
