function img=tz_synrgbimg(ch1,ch2,ch3,option)
%TZ_SYNRGBIMG Obsolete. See ML_SYNRGBIMG.
%   IMG = TZ_SYNRGBIMG(CH1,CH2,CH3,OPTION) returns a RGB image with three 
%   channels CH1,CH2 and CH3. One or two of the three channels could be
%   empty. But any two nonempty channels should have the same size.
%   There are three options for combination according to OPTION:
%       'org' - original values of the channels
%       'usc' - normalized over all channels
%       'isc' - normalized one by one

%   27-SEP-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

error(tz_genmsg('of','tz_synrgbimg','ml_synrgbimg'));

if isempty(ch1) & isempty(ch2) & isempty(ch3)
    img=[];
    warning('all channels are empty');
    return;
end

s1=size(ch1);
s2=size(ch2);
s3=size(ch3);

imgsize=max([s1;s2;s3]);

chs={ch1,ch2,ch3};

switch option
case 'org'
    for i=1:3
        if ~isempty(chs{i})
            img(:,:,i)=chs{i};
        else
            img(:,:,i)=zeros(imgsize);
        end
        
    end
case 'usc'
    minv=double(min(min([ch1;ch2;ch3])));
    maxv=double(max(max([ch1;ch2;ch3])));
    for i=1:3
        if ~isempty(chs{i})
            img(:,:,i)=tz_bcimg(double(chs{i}),[minv,maxv],[0 1]);
        else
            img(:,:,i)=zeros(imgsize);
        end    
    end
case 'isc'
    
    for i=1:3
        if ~isempty(chs{i})
            img(:,:,i)=tz_bcimg(double(chs{i}),[],[0 1]);
        else
            img(:,:,i)=zeros(imgsize);
        end    
    end
end
