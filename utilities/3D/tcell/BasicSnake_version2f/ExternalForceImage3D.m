function Eextern = ExternalForceImage3D(I,Wline, Wedge,Sigma, voxel_scales)
% Eextern = ExternalForceImage3D(I,Wline, Wedge,Sigma)
% 
% inputs, 
%  I : The image
%  Sigma : Sigma used to calculated image derivatives 
%  Wline : Attraction to lines, if negative to black lines otherwise white
%          lines
%  Wterm : Attraction to terminations of lines (end points) and corners
%
% outputs,
%  Eextern : The energy function described by the image
%
% Function is written by D.Kroon University of Twente (July 2010)

Ix=ImageDerivatives3D(I,Sigma,'x', voxel_scales);
Iy=ImageDerivatives3D(I,Sigma,'y', voxel_scales);
Iz=ImageDerivatives3D(I,Sigma,'z', voxel_scales);

Eline =  imgaussian(I,Sigma);
Eedge = sqrt(Ix.^2 + Iy.^2 + Iz.^2); 

% % 2013-01-14 tebuck: try to make lines recognized as laplacian responses:
% fprintf('>>>>>>>> HACK\n')
% Ixx=ImageDerivatives3D(I,Sigma,'xx');
% Iyy=ImageDerivatives3D(I,Sigma,'yy');
% Izz=ImageDerivatives3D(I,Sigma,'zz');
% Eline = -(Ixx + Iyy + Izz);
% % end mod


Eextern= (Wline*Eline - Wedge*Eedge); 

