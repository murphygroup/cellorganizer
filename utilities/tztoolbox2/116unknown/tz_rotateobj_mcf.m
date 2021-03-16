function newobjects = tz_rotatecombobj(dirname,combobjects,combclass,combcellidx)

%

load('/home/velliste/feat/class_names.mat');

NumberOfClasses = 10;
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = mv_classname( class);
    protfiles = mv_dir([dirname '/' class_name '/prot/*.dat']);
    cropfiles = mv_dir([dirname '/' class_name '/crop/*.tif*']);
    dnafiles = mv_dir([dirname '/' class_name '/dna/*.dat']);
    L = length( protfiles);
    if( (L ~= length( cropfiles)))
        error( 'There must be a matching number of protein images, DNA images and crop  images'); 
    end
    n_images = L;
    for N = 1 : n_images
        fprintf(1,[ class_name ', image ' num2str(N) '\n']);
        protimgname = [dirname '/' class_name '/prot/' protfiles{N}];
        cropimgname = [dirname '/' class_name '/crop/' cropfiles{N}];
        dnaimgname = [dirname '/' class_name '/dna/' dnafiles{N}];
        
        protimg=mv_readimage(protimgname);
        cropimg=mv_readimage(cropimgname);
        dnaimg=mv_readimage(dnaimgname);
        
        [procimg,prot_maskimage]=ml_preprocess(dnaimg,cropimg,'ml','yesbgsub');
        
        [theta,center]=tz_bwmajorangle(procimg);
        
        objs=combobjects(combclass==class & combcellidx==N);
        oldpos=tz_getobjcenter(objs);
        newpos=tz_rotate_2d([oldpos(:,1)-center(1),oldpos(:,2)-center(2)],-theta);
        
        drawnow
        newimg=imrotate(procimg,theta*180./pi);
        imshow(newimg,[]);
        mom=tz_bwmoment(newimg);
        hold on
        tz_plotgraph([newpos(:,1)+mom.cx,newpos(:,2)+mom.cy],[],1:2,'x','.',[]);
        hold off
        
        %cellobjs{N}=tz_findimgobjs(imagename,cropimagename,way);
    end
    %objects{class}=cellobjs;
end

% hold off
