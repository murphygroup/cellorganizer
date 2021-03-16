function objmodel=tz_bldobjmodel(cellobjcom,sel,ca)
%TZ_BLDOBJMODEL Build multinomial models from MCF object numbers.
%   OBJMODEL = TZ_BLDOBJMODEL(CELLOBJCOM,SEL)
%   
%   OBJMODEL = TZ_BLDOBJMODEL(CELLOBJCOM,SEL,CA)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function objmodel=tz_bldobjmodel(cellobjcom,sel,prenorm,ca)
%    
%OVERVIEW:
%   Build multinomial models from objnum MCF
%PARAMETERS:
%   cellobjcom - MCF for objnum
%   sel - select cells; useful for training and testing
%   ca - calibrate or not
%RETURN:
%   objmodel - multinomial model for each class, with average number per cell in the last column
%DESCRIPTION:
%   
%HISTORY:
%   03-MAR-2004 Initial write TINGZ
%   07-MAR-2004 Modified TINGZ
%   14-MAR-2004 Modified TINGZ
%   17-MAY-2004 Modified TINGZ
%       - add comments
%       - remove prenorm parameter

if ~exist('ca','var')
    ca=0;
end

nclass=length(cellobjcom);

for i=1:nclass
    if ~isempty(sel)
        objcom=cellobjcom{i}(sel,:);
    else
        objcom=cellobjcom{i};
    end
    
    avgobjnum=sum(objcom(:))/size(objcom,1);

    objmodel(i,:)=[tz_normobjcom(sum(objcom,1),ca),avgobjnum];
end
