function normobjcom=tz_normobjcom(cellobjcom,ca,label)

% normobjcom=tz_normobjcom(cellobjcom,ca,label)
%
%OVERVIEW:
%   Normalize object compositions to multinomial distributions
%PARAMETERS:
%   cellobjcom - mxn Matrix for n object compositions for m cells
%   ca - calibration flag. to avoid zero probabilities
%   label - tell if the parameter cellobjcom has labels in its last column
%RETURN:
%   normobjcom - multinomial parameters
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   18-APR-2004 Modfied TINGZ
%       - Add calibratioin
%   11-Mar-2004 Modfied TINGZ
%       - Remove unlabel bugs
%   05-NOV-2004 Modfied TINGZ
%       - revise comments

if ~exist('label','var')
    label=0;
end

if ~exist('ca','var')
    ca=0;
end

if label
    ulcellobjcom=cellobjcom(:,1:(end-1));
else
    ulcellobjcom=cellobjcom;
end

if ca
   %Calibration 
   %The probability of 0 is 0.95
   ncell=size(ulcellobjcom,1);
   N=sum(ulcellobjcom,2);
   cavalue=1-0.95.^(1./N);
   for i=1:ncell
       ulcellobjcom(i,ulcellobjcom(i,:)==0)=cavalue(i);
   end
end

if label==0
    normobjcom=ulcellobjcom./repmat(sum(ulcellobjcom,2),1,size(ulcellobjcom,2));
else
    normobjcom=[ulcellobjcom./repmat(sum(ulcellobjcom,2),1,size(ulcellobjcom,2)),cellobjcom(:,end)];
end


