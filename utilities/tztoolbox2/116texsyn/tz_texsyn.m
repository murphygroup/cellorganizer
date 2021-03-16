function img2=tz_texsyn(img,mask,param)
%TZ_TEXSYN Synthesize texture.
%   IMG2 = TZ_TEXSYN(IMG,MASK) returns an image that is the synthesized
%   texture similar with the texture of IMG. MASK a the binary image for
%   masking background. If there is no mask, set MASK to [].
%   
%   IMG2 = TZ_TEXSYN(IMG,MASK,PARAM) allows selecting texture synthesis
%   method and specifying parameters. PARAM is a structure with the
%   following fields:
%       'method' - synthesis method
%           'nb1' - nearest neighbor (too slow)
%           'com' - combination of nearest neighbor and hybrid texture
%               synthesis.
%           'hyb' - hybrid texture synthesis.
%           'wav' - pyramid texture synthesis
%
%   See also

%   12-Jan-2006 Initial write T. Zhao
%   14-May-2006 Modeified T. Zhao
%       - add 'hyb' option. Change old 'hyb' to 'com'.
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('method','hyb','texsize',[256,256]));

methodsUsingHyb = {'hyb','com'};

if ismember(param.method,methodsUsingHyb)
    param = ml_initparam(param,struct('cropbg',0,'errtol',0.5, ...
        'patchsize',[32 32]));

    if ~isempty(mask)
        img = tz_imcropbg2(img,mask);
        mask = tz_imcropbg2(mask,mask);
    else
        if param.cropbg==1
            img = ml_imcropbg(img);
            mask = img>0;
        end
    end

    dimg=img/max(img(:));

    data = ...
        tz_gen_input(dimg,param.texsize(1),param.texsize(2), ...
        param.errtol,1,'nowrap',1.0,1.0,...
        param.patchsize(1),param.patchsize(2),'simple',mask);
    [knn, distance] = build_knn(dimg,7,7,97,1);
    out=hybridsynthesize(data,7,knn);
end

switch param.method
    case 'nb1'
        param = ml_initparam(param, ...
            struct('k',5,'ep',5,'df','pixel'));
        
        k = param.k;
        ep = param.ep;
        df = param.df
        
        initpixels=tz_genimgpixel(img,mask,(k+1).^2);

        img2mask=zeros(size(img));
        img2=img2mask;

        img2(1:k+1,1:k+1)=reshape(initpixels,k+1,k+1);
        img2mask(1:k+1,1:k+1)=1;
        img2mask(k+1,k+1)=0;

        for i=1:size(img,1)
            for j=1:size(img,2)
                if img2mask(i,j)==0
                    wnd=tz_getnbwnd(img2,[i,j],k,1);
                    qmask=tz_getnbwnd(img2mask,[i,j],k,1);
                    rpixel=tz_genpixel(img,mask,wnd,qmask,df,Inf,ep)
                    img2mask(i,j)=1;
                    img2(i,j)=rpixel;
                end
            end
        end
    case 'hyb'
        img2 = out(:,:,1);
        img2 = img2.*max(img(:));
    case 'com'
        param = ml_initparam(param,struct('Nsc',4,'Nor',4,'Na',7,'Niter',50));
        
        Nsc = param.Nsc;
        Nor = param.Nor;
        Na = param.Na; 

        params = textureAnalysis(out(:,:,1), Nsc, Nor, Na);

        Niter = param.Niter;	
%         Nsx = param.texsize(1);	
%         Nsy = param.texsize(2);	

        figure

        [img2,snrP,imS]= textureSynthesis(params, param.texsize, Niter);
    case 'wav' 
        param = ml_initparam(param,struct('Nsc',4,'Nor',4,'Na',7,'Niter',50));
        
        Nsc = param.Nsc; % Number of scales
        Nor = param.Nor; % Number of orientations
        Na = param.Na;  % Spatial neighborhood, must be an odd number

        params = textureAnalysis(img, Nsc, Nor, Na);

        Niter = param.Niter;	% Number of iterations of synthesis loop
        Nsx = param.texsize(1);	% Size of synthetic image is Nsy x Nsx
        Nsy = param.texsize(2);	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

        figure

        [img2,snrP,imS]= textureSynthesis(params, [Nsy Nsx], Niter);
        
end

