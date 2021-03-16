% Example 1: Synthesis of a "text" texture image, using
% Portilla-Simoncelli texture analysis/synthesis code, based on
% alternate projections onto statistical constraints in a complex
% overcomplete wavelet representation.
%
% See Readme.txt, and headers of textureAnalysis.m and
% textureSynthesis.m for more details.
%
% Javier Portilla (javier@decsai.ugr.es).  March, 2001

% [bigobj,classidx]=tz_findbigobj(combobjects,combclass,1000);

%close all

%im0 = pgmRead('metal.pgm');	% im0 is a double float matrix!
% im0=imread('/home/tingz/matlab/data/mhela/erdak/mhelacell_erdak10.png');
% %im0=im0(127:382,257:512);
% im0=double(im0(1:1024,1:1024));
%sel=1;
%im0 = tz_obj2rect(bigobj{sel},[256 256]); %
%im0=im0(1:30,1:30);
%img=tz_objimg(combobjects{201},[]);
% im0=img(30:60,30:60);

%img=imread('');
%im0=double(proj_image(30:200,30:200));

selobj=bigobj(classidx==10);
im0=tz_objimg(selobj{1},[]);
% %im0=im0(1:128,1:128);
im0=imresize(im0,[256 256],'bilinear');
%im0=img;
%im0=im0(1:64,1:64);

Nsc = 4; % Number of scales
Nor = 4; % Number of orientations
Na = 3;  % Spatial neighborhood is Na x Na coefficients
	 % It must be an odd number!

params = textureAnalysis(im0, Nsc, Nor, Na);

Niter = 50;	% Number of iterations of synthesis loop
Nsx = 256;	% Size of synthetic image is Nsy x Nsx
Nsy = 256;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

[res,snrP,imS]= textureSynthesis(params, [Nsy Nsx], Niter);

close all
figure(1)
showIm(im0, 'auto', 1, 'Original texture');
figure(2)
showIm(res, 'auto', 1, 'Synthesized texture 1');
figure(3)
showIm(tz_reorgimg(res),'auto',1,'Synthesized texture');

% Can you read the NEW text? ;-)

% figure
% img=tz_objimg(bigobj{sel},[]);
% img=imresize(img,[256 256],'bilinear');
% showIm(img,'auto',1,'original object')
% 
% figure
% showIm(tz_rect2obj(tz_reorgimg(im0),bigobj{sel}),'auto',1,'original object 2')
% figure
% showIm(tz_rect2obj(tz_reorgimg(res),bigobj{sel}),'auto',1,'reconstructed object')
