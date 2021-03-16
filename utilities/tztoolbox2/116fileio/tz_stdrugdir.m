function tz_stdrugdir(datadir)
%TZ_STDRUGDIR Standarize drug directory for batch processing.
%   TZ_STDRUGDIR(DATADIR) is used before batch processing drug
%   images acquired from QED. DATADIR is where images are saved.

%   ??-???-2004 Initial write T. Zhao
%   03-NOV-2004 Modified T. Zhao
%       - add comments
%   19-NOV-2004 Modified T. Zhao
%       - automatic dir finding
%   Copyright (c) Murphy Lab, Carnegie Mellon University

warning('dir changed');

protlist=ml_dir(datadir);
protlist=tz_cleandirs(protlist);
nprot=length(protlist);

for i=1:nprot
    protdir=[datadir '/' protlist{i}];
    druglist=ml_dir(protdir);
    druglist=tz_cleandirs(druglist);
    ndrug=length(druglist);
    
    for j=1:ndrug
        drugdir=[protdir '/' druglist{j}];
        celllist=ml_dir(drugdir);
        celllist=tz_cleandirs(celllist);
        ncell=length(celllist);
        
        for k=1:ncell
            celldir=[drugdir '/' celllist{k}];
            tmpdirlist=ml_dir(celldir);
            tmpdirlist=tz_cleandirs(tmpdirlist);
            for i=1:length(tmpdirlist)
                if all(tmpdirlist{i}(end-3:end)=='Data')
                    sampledir=[celldir '/' tmpdirlist{i}];
                    %nsampledir=[celldir '/' '~sample Data'];
                    if exist(sampledir,'dir')
                        spos=find(sampledir==' ');
                        sampledir=tz_insertvec(sampledir,'\',spos);
                        %sampledir=[celldir '/' 'Movie\ Data'];
                        unix(['mv ' sampledir '/*.tif ' celldir]);
                    end
                end
            end
            
%             if exist(nsampledir,'dir')
%                 sampledir=[celldir '/' '~sample\ Data'];
%                 unix(['mv ' sampledir '/*.tif ' celldir]);
%             end
            
        end
        unix(['rm ' drugdir '/.DS*']);
    end
end
