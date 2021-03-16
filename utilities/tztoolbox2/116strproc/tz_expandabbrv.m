function [fullword,stat] = tz_expandabbrv(s)
%TZ_EXPANDABBRV Expansion of an abbreviation.
%   TZ_EXPANDABBRV(S) returns the full word for which S stands. It is
%   useful for unstanding the codes in this toolbox.

%   03-Aug-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

stat = 1;

switch (s)
    case '2'
        fullword = 'to';
    case 'abbrv'
        fullword='abbrevation';
    case 'arg'
        fullword='argument';
    case 'classif'
        fullword='classification';
    case 'clst'
        fullword='clustering';
    case 'cmt'
        fullword='comment';
    case 'comb'
        fullword='combine';
    case 'comp'
        fullword='comparison';
    case 'datafeat'
        fullword='data features';
    case 'dataproc'
        fullword='data processing';
    case 'datastr'
        fullword='data structure';
    case 'dep'
        fullword='dependent';
    case 'eval'
        fullword='evaluate';
    case 'feat'
        fullword='feature';
    case 'fila'
        fullword='filament';
    case 'fileio'
        fullword='file I/O';
    case 'fun'
        fullword='function';
    case 'gen'
        fullword='generate';
    case 'geom'
        fullword='geometry';
    case 'graph'
        fullword='graphics';
    case 'img'
        fullword='image';
    case 'imgcode'
        fullword='image coding';
    case 'imgfeat'
        fullword='image features';
    case 'imgproc'
        fullword='image processing';
    case 'imgstat'
        fullword='image statistics';
    case 'imgsyn'
        fullword='image synthesis';
    case 'init'
        fullword='initialize';
    case 'math'
        fullword='basic mathematics';
    case 'msg'
        fullword='message';
    case 'mt'
        fullword='matlab';
    case 'num'
        fullword='number';
    case 'obj'
        fullword='object';
    case 'objdetect'
        fullword='object detection';
    case 'objfeat'
        fullword='object features';
    case 'objproc'
        fullword='object processing';
    case 'parallel'
        fullword='parallel processing';
    case 'pos'
        fullword='position';
    case 'pt'
        fullword='point';
    case 'reg'
        fullword='regress';
    case 'stat'
        fullword='statistics';
    case 'str'
        fullword='string';
    case 'strproc'
        fullword='string processing';
    case 'texsyn'
        fullword='texture synthesis';
    case 'tok'
        fullword='token';
    case 'tz'
        fullword='Ting Zhao';
    case 'mxs'
        fullword = 'medial axis';
    case 'crd'
        fullword = 'coordinate';
    case 'vis3d'
        fullword = '3D visualization';
    otherwise
        fullword='unknown word';
        stat = 0;
end
