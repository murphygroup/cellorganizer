function ras=tz_rarea_mcf(imgdir)
%TZ_RAREA_MCF Calcualte the area ratios between cells and their nuclei.
%   RAS = TZ_RAREA_MCF(IMGDIR) returns the vector of area ratios between
%   cells and their nuclei. The cells are loaded from the directory IMGDIR.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function tz_parsecell_mcf(imgdir,savedir)
%

combcodes={};
combclass=[];
combcellidx=[];

files=dir(imgdir);
classes={};
for i=1:length(files)
    if isdir([imgdir '/' files(i).name]) & files(i).name(1)~='.'
        classes={classes{:},files(i).name};
    end
end
channels={'prot','dna','cell'};
imgsize=[1024,1024];
ras=[];
for pi=1:length(classes)
    patterndir=[imgdir '/' classes{pi}];
    protdir=[patterndir '/' channels{1}];
   
    cellbodydir=[patterndir '/cellbody'];
    nucbodydir=[patterndir '/nucbody'];
    prots=dir([protdir '/*.mat']);
    

    for ci=1:length(prots)
        [pi ci]
        [tmp,dpos]=find(prots(ci).name=='_');
        cellname=prots(ci).name(1:dpos-1);
        
        load([cellbodydir '/' cellname '_cell.mat']);
        load([nucbodydir '/' cellname '_dna.mat']);         
        ras=[ras,length(cellbody(:))/length(nucbody(:))];
    end
end