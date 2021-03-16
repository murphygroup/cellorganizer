function Energy = AverageEnergy(segcell)
% grj 10/22/14 - set contour visibility to 'off'
% warning off

% folder = './FAAS_sample/ThreshCellMasks/';
% folder = '../detailed_cell_shape_interpolation_for_devin/2Dmasks/';
% folder = '../CellShape/Cyto116b/ts040629a001/';

% num = length(files);
%ask gong, why is this set to 8? 
KNN_half = 8;
%This should also probably not be hard coded, but is only being used for
%the histogram right now so should be ok. 
num_bin = 256;
% result = cell(num,1);
st=zeros(1,num_bin);
value=st;

obj = squeeze(sum(segcell,3))>0;

co = contour(obj,[0.5,0.5], 'Visible', 'off');


co = co';
co = co(2:end,:);
N = size(co,1);

curvature = zeros(N,1);


for j=1:N
    curvature(j) = est_curv(obj,co,j,N,KNN_half);
end
%     figure(2),imshow(zeros(size(obj)));hold on
%     co2 = co';
%     surface([co2(1,:);co2(1,:)],[co2(2,:);co2(2,:)],[curvature';curvature'],'facecol','no','linew',10)

Energy = sum(curvature.^2);
% [st(i,:),value(i,:)]=hist(curvature{i},num_bin);
% result{i}.st = st(i,:);
% result{i}.value = value(i,:);


