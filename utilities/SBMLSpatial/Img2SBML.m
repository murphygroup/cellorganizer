function Img2SBML(imgpath,names,ordinalList,savepath,resolution,SBMLfile)
%This function takes in an image path and a name creates an SBML-spatial
%file for it. 


if ischar(imgpath)
    imgpath = {imgpath};
end

if ischar(names)
    names = {names};
end

if nargin<3
    warning('No ordinals specified. Assuming all ordinals are 0');
    ordinalList = [];
end
if nargin<4
    savepath = './model.xml';
end
if nargin<5
    warning('No resolution specified. Not rescaling image.');
    resolution = [1,1,1];
end
if nargin<6
    SBMLfile = [];
end



if length(imgpath)~=length(names)
    error('Number of images does not match number of compartment names!')
end

CSGboxes = struct;
CSGboxes.class = ['BoundingBoxes'];
CSGboxes.list = [];
%loop through each image, read it in and make the SBML-spatial file. 
nc = 1;%keep track of the number of compartments (nc)
for i = 1:length(names)
    %read in the image 
    img = ml_readimage(imgpath{i});
    for z = 1:size(img,3)
        img(:,:,z) = imfill(img(:,:,z),'hole');
    end
    levels = unique(img);
    %if the image is not a binary image, try to find names for each level
    if length(levels)==1
        error('Only one pixel level present in image!');
    elseif length(levels)>2
        for j = 2:length(levels)
            %assign the compartment name.
            param.SBML_Name{nc} = names{i}{j};
            imgs{nc} = img==levels(j);
            param.ind = nc;
            CSGboxes = getBox(imgs{nc},['BB',names{i}],resolution,CSGboxes);
%             CSGboxes.list(nc) = getBox(imgs{nc},names{i},[1,1,1],param);
            if isempty(ordinalList)
                param.ordinal{nc} = 0;
            else
                param.ordinal{nc} = ordinalList{i}{j};
            end
            nc = nc+1;
        end
        
    else
        imgs{nc} = img>levels(1);
        param.SBML_Name{nc} = names{i};
        param.ind = nc;
        CSGboxes = getBox(imgs{nc},['BB',names{i}],resolution,CSGboxes);
%         CSGboxes.list(nc) = getBox(imgs{nc},names{i},[1,1,1],param);
        if isempty(ordinalList)
            param.ordinal{nc} = 0;
        else
            param.ordinal{nc} = ordinalList{i};
        end
        nc = nc+1;
        

    end
    
end


sbmlStruct = createImg2SBMLstruct(imgs,param);
% CSGboxes.primitiveOnly = 1;

[ result ] = instance2SBML( CSGboxes, sbmlStruct, savepath, SBMLfile,resolution)


end 

