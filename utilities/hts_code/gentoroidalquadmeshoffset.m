function [V,F,BB,rows,cols,VF] =...
    gentorodialquadmeshoffset(dy,dx,numrows,numcols,rowoff,coloff)
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
% File gentoroidalquadmeshoffset.m
%   This subroutine will compute a uniform quad mesh
%   with mesh patches of size dy x dx and numx,numy patches
%   in x,y-dimension. the mesh, although in integer 2D
%   coordinates, wraps in 2 dimensions (toroidal topology),
%   and the wrapping patches are explicitly given by the
%   bounding box coordinates
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
%   VF - [vertex -> face] adjacency list (nv x 4)
%
%   the information returned can be used to traverse the
%   triangle mesh along the ordering within the face list
%

% all loops below can and should be merged, but are currently kept separate for
% readability

% start by building the vertex list
V = zeros(numrows*numcols,2);
pos = 0;
for j=1:numrows,
    for i=1:numcols,
        pos = pos + 1;           % current position in the vertex list
        V(pos,1) = rowoff + (j-1)*dy + 1; % row pixel coord
        V(pos,2) = coloff + (i-1)*dx + 1; % column pixel coord
        % make wrap
        if (V(pos,1) > (numrows*dy)), V(pos,1) = V(pos,1) - (numrows*dy); end
        if (V(pos,2) > (numcols*dx)), V(pos,2) = V(pos,2) - (numcols*dx); end
    end
end

% now build the face list F and bounding box BB information
F  = zeros(numrows*numcols,4);
% we store BB information explicitly in case we generalize the 
% polygonal patches, which will be one goal of this project
BB = zeros(numrows*numcols,4);

pos = 0;
for j=1:numrows,
    for i=1:numcols,
        % current position in the face list F
        pos = pos + 1;

        % upper left face corner...
        F(pos,1) = pos; 
        % ...which is also BBmin
        BB(pos,1) = V(F(pos,1),1); 
        BB(pos,2) = V(F(pos,1),2);

        % upper right corner
        if rem(pos,numcols) == 0, 
            % pos = last position in this row, wrap it (make it first position
            % in the current row)
            F(pos,2) = (j-1)*numcols + 1;
        else
            % next position to the right
            F(pos,2) = pos + 1;
        end

        % lower right corner (start from upper right corner = F(pos,2))...
        if (F(pos,2) + numcols) > (numrows*numcols), 
            % pos = last position in this column, wrap it (make it first position
            % in the current column)
            if (i + 1 > numcols), F(pos,3) = rem(i + 1, numcols); else
                F(pos,3) = i + 1;  end
        else
            % next position below
            F(pos,3) = F(pos,2) + numcols;
        end
        % ... which is also BBmax
        BB(pos,3) = V(F(pos,3),1); 
        BB(pos,4) = V(F(pos,3),2);
        
        % lower left corner
        if (pos + numcols) > (numrows*numcols),
            % pos = last position in this column, wrap it (make it first position
            % in the current column)
            F(pos,4) = i;
        else
            % next position below
            F(pos,4) = pos + numcols;
        end

    end
end

rows = dy * numrows;
cols = dx * numcols;

% build VF list. optionally used for recomputing the BB coords when 
% perturbing the vertex coordinates
VF = zeros(numrows*numcols,2);
pos = 0;
for j=1:numrows,
    for i=1:numcols,
        % current position in the vertex list V
        pos = pos + 1;
        
        % bottom right face
        VF(pos,1) = pos;
        
        % bottom left face
        if (rem(pos-1,numcols) == 0), VF(pos,2) = pos-1+numcols; 
        else VF(pos,2) = pos - 1; end
        
        % upper right face
        if ((pos - numcols) < 1), VF(pos,3) = numrows*numcols+(pos-numcols);
        else VF(pos,3) = pos - numcols; end
        
        % upper left face (start from upper right face)
        if (rem(VF(pos,3) - 1,numcols) == 0), VF(pos,4) = VF(pos,3) - 1 + numcols;
        else VF(pos,4) = VF(pos,3) - 1; end
    end
end
