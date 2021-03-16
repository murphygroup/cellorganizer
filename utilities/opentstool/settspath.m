function settspath(TSTOOLpath)

TSTOOLpath=fullfile(TSTOOLpath,'tstoolbox');
path(TSTOOLpath, path);
path(fullfile(TSTOOLpath, 'demos'), path);
path(fullfile(TSTOOLpath, 'gui'), path);
path(fullfile(TSTOOLpath, 'utils'), path);
path(fullfile(TSTOOLpath, 'mex'), path);
path(fullfile(TSTOOLpath, fullfile('mex', mexext)), path);

