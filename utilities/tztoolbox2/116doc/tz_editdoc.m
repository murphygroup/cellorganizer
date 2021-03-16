function tz_editdoc(category,docname)
%TZ_EDITDOC Edit documents.
%   TZ_EDITDOC(CATEGORY,DOCNAME) provides a shortcut for opening a document
%   to edit. CATEGORY is the category of the target document:
%       'main' - main categories. The corresponding DOCNAME could be
%           'index' - doc/index.html
%           'log' - doc/log.html
%           'glossary' - doc/glossary.html
%           'guide' - doc/userguide/ug_index.html
%           'start' - doc/learn/gs_index.html
%       'proj' - projects. The corresponding DOCNAME could be
%           'index' - doc/projects/pj_index.html
%           'genmodel' - doc/projects/genmodel.html
%           'yeast' - doc/projects/yeast.html
%       'dep' - developemnt. The corresponding DOCNAME could be
%           'index' - doc/development/dp_index.html
%           'design' - doc/development/design.html
%           'file' - doc/development/filestruct.html
%           'cmt' -  writing comments
%           'imgproc' - image processing
%           'stat' - statistics
%           'math' - mathematics
%           'model' - modeling
%   
%   See also

%   16-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

global glMtpath glShareddir glDocdir

docRootDir = [glMtpath filesep glShareddir filesep glDocdir];
editor = 'emacs';

switch category
    case 'main'
        switch docname
            case 'index'
                targetFile = 'index.html';
            case 'log'
                targetFile = 'log.html';
            case 'glossary' 
                targetFile = 'glossary.html';
            case 'guide'
                targetFile = 'userguide/ug_index.html';
            case 'start'
                targetFile = 'learn/gs_index.html';
            case 'test'
                targetFile = 'testing/index.html';
            otherwise
                error(['Unrecognized doc name in the category main. ' ...
                    'It must be one of the following options: ' ...
                    '''index'', ''log'', ''glossary'', ''guide'',' ...
                    '''start'',''test''']);
        end
    case 'proj'
        switch docname
            case 'index' 
                targetFile = 'projects/pj_index.html';
            case 'genmodel'
                targetFile = 'projects/genmodel.html';
            case 'yeast'
                targetFile = 'projects/yeast.html';
            otherwise
                error('Unrecognized doc name in the category project');
        end
    case 'dep'
        switch docname
            case 'index' 
                targetFile = 'development/dp_index.html';
            case 'design' 
                targetFile = 'development/design.html';
            case 'file'
                targetFile = 'development/filestruct.html';
            case 'cmt'
                targetFile = 'development/mtcmt.html';
            case 'imgproc'
                targetFile = 'development/imageproc.html';
            case 'stat'
                targetFile = 'development/stat.html';
            case 'math'
                targetFile = 'development/math.html';
            case 'model'
                targetFile = 'development/model.html';
            otherwise
                error(['Unrecognized doc name in the category development.' ...
                       'It must be ''index'' or ''design''.']);           
        end
    otherwise
        error('Unrecognized category.'); 
end

targetPath = [docRootDir filesep targetFile];
disp(['Opening ' targetPath '...']);
system([editor ' ' targetPath ' &']);
