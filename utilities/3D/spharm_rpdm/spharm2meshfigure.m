function [Zvert, fs] = spharm2meshfigure(deg,fvec,meshtype,plot,figtitle,filename,dpi)
% generate 3D mesh and optional figure from a spherical parameterization
% R.F.Murphy 09/12/2022
% Modified 9/14/2022 to allow specification of mesh resolutions in meshtype

% inputs
% deg = degree of spharm descriptors
% fvec = spharm descriptors
% meshtype.type = (optional) 'even' or 'triangular' (default 'even')
% meshtype.nPhi = (optional) number of meridians (if type='even')
% meshtype.nTheta = (optional) number of parallels (if type='even')
% meshtype.nVerticies = (optional) number of vertices (if type='triangular')
% plot = (optional) generate figure if true (default false)
% figtitle = (optional) title to put on figure
% filename = (optional) filenname to save figure to (default not saved)
% dpi = (optional) dots per inch % figtitle = (optional) title to put on figure (default none)
%
% output
% Zvert = vertices
% fs = faces
%

    if ~exist('meshtype','var')
        meshtype = [];
        meshtype.type = 'even';
    end
    if strcmp(meshtype.type,'even')
        % generate evenly distributed vertices (vs) on the surface of a 
        % sphere and their connectedness (faces, fs)
        if ~isfield(meshtype,'nPhi') meshtype.nPhi = []; end
        if ~isfield(meshtype,'nTheta') meshtype.nTheta = []; end
        [vs fs] = sphereMesh([0 0 0 1 meshtype.nPhi meshtype.nTheta]);
    else
        % generate a triangular mesh
        if ~isfield(meshtype,'nVertices') meshtype.nVertices = 4002; end
        [vs, fs]=SpiralSampleSphere(meshtype.nVertices);
    end
    % convert vertices to SPHARM basis function values
    Zs = calculate_SPHARM_basis(vs, deg);
    % adjust the basis functions using the spherical parameterization
    Zvert = real(Zs*fvec);
    if exist('plot', 'var') && plot
        if ~exist('figtitle','var') figtitle = []; end
        if ~exist('filename','var') filename = []; end
        if ~exist('dpi','var') dpi = []; end
        mesh2figure(Zvert,fs,figtitle,filename,dpi)
    end
    close
end
