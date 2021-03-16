function [nucpts,cellpts] = ...
    tz_synucell(cellpts,nucpts,sizeratio,angle,move)
%TZ_SYNUCELL Synthesize a cell containing a nucleus.
%   NUCPTS = TZ_SYNUCELL(CELLPTS,NUCPTS,SIZERATIO,ANGLE,MOVE) returns the
%   boundary points of nucleus according to the size ratio SIZERATIO, 
%   relative angle ANGLE and offset MOVE between the cell CELLPTS and
%   nucleus NUCPTS.
%   
%   [NUCPTS,CELLPTS] = TZ_SYNUCELL(...) also returns the boundary points of
%   the cell.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function pts = tz_synucell(cellpts,nucpts,angle,move)
%OVERVIEW
%   
%PARAMETERS
%   cellpts - 
%   nucpts - 
%   angle - 
%   move - 
%RETURN
%   pts - 
%DESCRIPTION
%   
%HISTORY
%   25-Apr-2005 Initial write TINGZ
%SEE ALSO
%   


[cellcenter,cellmangle] = tz_edgecenter(cellpts);
[nuccenter,nucmangle] = tz_edgecenter(nucpts);

theta=cellmangle;
majorvec=[cos(theta),sin(theta)];
mproj=cellpts*majorvec';
cellmajor=max(mproj)-min(mproj);

minorvec=[sin(theta),-cos(theta)];
mproj=cellpts*minorvec';
cellminor=max(mproj)-min(mproj);
    
nuctheta=nucmangle;
nucmajorvec=[cos(nuctheta),sin(nuctheta)];
nucmproj=nucpts*nucmajorvec';
nucmajor=max(nucmproj)-min(nucmproj);

cellpts=round(cellpts*nucmajor*sizeratio/cellmajor);
[cellcenter,cellmangle] = tz_edgecenter(cellpts);
theta=cellmangle;
majorvec=[cos(theta),sin(theta)];
mproj=cellpts*majorvec';
cellmajor=max(mproj)-min(mproj);

minorvec=[sin(theta),-cos(theta)];
mproj=cellpts*minorvec';
cellminor=max(mproj)-min(mproj);
radiff=cellmangle-nucmangle+angle*pi/180;

pts=tz_rotate_2d(nucpts,radiff);
pts=ml_addrow(pts,-mean(pts,1));

tnuccenter=(move.*[cellmajor cellminor])*[majorvec;minorvec];
tnuccenter=cellcenter+tnuccenter;
nucpts=ml_addrow(pts,tnuccenter);