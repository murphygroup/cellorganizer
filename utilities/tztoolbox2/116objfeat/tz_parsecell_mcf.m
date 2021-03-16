function [combcodes,combclass,combcellidx]=tz_parsecell_mcf(imgdir,savedir)
%TZ_PARSECELL_MCF Obsolete. See TZ_PARSECELL2_MCF.
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
for pi=1:length(classes)
    
    patterndir=[imgdir '/' classes{pi}];
    protdir=[patterndir '/' channels{1}];
   
    cellbodydir=[patterndir '/cellbody'];
    nucbodydir=[patterndir '/nucbody'];
    prots=dir([protdir '/*.mat']);
    
    fullsavedir=[patterndir '/' savedir];
    
    if ~exist(fullsavedir,'dir')
        mkdir(patterndir,savedir);
    end
        
    for ci=1:length(prots)
        [pi ci]
        [tmp,dpos]=find(prots(ci).name=='_');
        cellname=prots(ci).name(1:dpos-1);
        savefile=[fullsavedir '/' cellname '.mat'];
        
        if ~exist(savefile,'file')            
            load([cellbodydir '/' cellname '_cell.mat']);
            load([nucbodydir '/' cellname '_dna.mat']);         
            
            [cellcontour,nuccontour,celldist,nucdist,cellhitpts,nuchitpts,angles,nuccenter]=...
                tz_parsecell(cellbody,nucbody,1,imgsize,0);
            save(savefile,'cellcontour','nuccontour',...
                'celldist','nucdist','cellhitpts','nuchitpts','angles','nuccenter');
        else
            load(savefile);
        end
        codes=struct('cellcontour',cellcontour,'nuccontour',nuccontour,...
                'celldist',celldist,'nucdist',nucdist,'cellhitpts',cellhitpts,...
                'nuchitpts',nuchitpts,'angles',angles,'nuccenter',nuccenter);
        combcodes={combcodes{:},codes};
        combclass=[combclass;pi];
        combcellidx=[combcellidx;ci];
    end
end