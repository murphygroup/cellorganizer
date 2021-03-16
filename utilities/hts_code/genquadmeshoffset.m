function [V,F,BB,rows,cols] =...
    genquadmeshoffset(dy,dx,numrows,numcols,rowoff,coloff)
%
% |----------------------------------------------------------|
% | Hybrid Texture Synthesis MATLAB package                  |
% |                                                          |
% | Author: Andrew Nealen                                    |
% |         Discrete Geometric Modeling Group                |
% |         Technische Universitaet Darmstadt, Germany       |
% |                                                          |
% | Note:   This is part of the prototype implementation     |
% |         accompanying our paper/my thesis                 |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen and M. Alexa |
% |         Eurographics Symposium on Rendering 2003         |
% |                                                          |
% |         Hybrid Texture Synthesis. A. Nealen              |
% |         Diplomarbeit (MSc), TU Darmstadt, 2003           |
% |                                                          |
% |         See the paper/thesis for further details.        |
% |----------------------------------------------------------|
%
% File genquadmeshoffset.m
%   This subroutine will compute a uniform quad mesh
%   with mesh patches of size dy x dx and numx,numy patches
%   in x,y-dimension. 
%
%   Input Values: (y = rows, x = cols)
%   dy, dx - width of each mesh patch (in pixels) in y and x dimension
%   numrows, numcols - number of mesh patches in y and x dimension
%   rowoff, coloff - offset in row and col pixels
%
%   Output Values:
%   V  - the vertex list (dim = nv x 2), 2D integer coordinates
%   F  - the face list   (dim = nf x 4), indices into vertex list.
%        starts upper left and runs clockwise around each face
%   BB - bounding box of each patch (dim = nf x 4)
%        the index corresponds to the index into the facelist
%        with:
%        BB(i,1) = BBmin_y, BB(i,2) = BBmin_x, and
%        BB(i,3) = BBmax_y, BB(i,4) = BBmax_x
%        the bounding box ccords are considered 'points of
%        reference' in the texture synthesis algorithm
%        (where the patch may be deformed during synthesis)
%   rows - number of pixel rows in synthesis result
%   cols - number of pixel cols in synthesis result
%
%   the information returned can be used to traverse the
%   triangle mesh along the ordering within the face list
%

% all loops below can and should be merged, but are currently kept separate for
% readability

% start by building the vertex list
V = zeros((numrows+1)*(numcols+1),2);
pos = 0;
for j=1:(numrows+1),
    for i=1:(numcols+1),
        pos = pos + 1;           % current position in the vertex list
        V(pos,1) = rowoff + (j-1)*dy + 1; % row pixel coord
        V(pos,2) = coloff + (i-1)*dx + 1; % column pixel coord
    end
end

% now build the face list F and bounding box BB information
F  = zeros(numrows*numcols,4);
% we store BB information explicitly. for (possible) generalization purposes
BB = zeros(numrows*numcols,4);

pos = 0;
for j=1:numrows,
    for i=1:numcols,
        % current position in the face list F
        pos = pos + 1;
        rowindex = ceil(pos/numcols);
        colindex = pos - (rowindex-1)*numcols;
        
        % upper left face corner...
        F(pos,1) = (rowindex-1)*(numcols+1)+colindex; 
        % ...which is also BBmin
        BB(pos,1) = V(F(pos,1),1); 
        BB(pos,2) = V(F(pos,1),2);

        % upper right corner, next position to the right
        F(pos,2) = F(pos,1) + 1;

        % lower right corner (start from upper right corner = F(pos,2))...
        % next position below
        F(pos,3) = F(pos,2) + (numcols+1);

        % ... which is also BBmax
        BB(pos,3) = V(F(pos,3),1); 
        BB(pos,4) = V(F(pos,3),2);
        
        % lower left corner, next position below
        F(pos,4) = F(pos,1) + (numcols+1);
    end
end

rows = dy * numrows;
cols = dx * numcols;
