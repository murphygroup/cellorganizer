function distribution_tif(filepath,savepath)
warning off

% folder = './FAAS_sample/ThreshCellMasks/';
% folder = '../detailed_cell_shape_interpolation_for_devin/2Dmasks/';
% folder = '../CellShape/Cyto116b/ts040629a001/';
files = ml_dir([filepath filesep '*.tif']);
num = length(files);
%ask gong, why is this set to 8? 
KNN_half = 8;
%This should also probably not be hard coded, but is only being used for
%the histogram right now so should be ok. 
num_bin = 256;
result = cell(num,1);
st=zeros(num,num_bin);
value=st;
Energy = zeros(num,1);
figure(10)
set(gca,'ColorOrder',jet(num))
colors = get(gca,'colororder');
curvature = cell(1,num);
for i=1:num

%     load(strcat(folder,'cell',num2str(i)));
    currentfile = [filepath filesep files{i}]
%     load(currentfile);
    segcell = ml_readimage(currentfile);
    
    obj = squeeze(sum(segcell,3));

    figure(1),co = contour(obj,[0.5,0.5],'color',colors(i,:));hold on
    
    
    co = co';
    co = co(2:end,:);
    N = size(co,1);

    curvature{i} = zeros(N,1);
    

    for j=1:N
        curvature{i}(j) = est_curv(obj,co,j,N,KNN_half);
    end
%     figure(2),imshow(zeros(size(obj)));hold on
%     co2 = co';
%     surface([co2(1,:);co2(1,:)],[co2(2,:);co2(2,:)],[curvature';curvature'],'facecol','no','linew',10)

    Energy(i) = sum(curvature{i}.^2);
    [st(i,:),value(i,:)]=hist(curvature{i},num_bin);
    result{i}.st = st(i,:);
    result{i}.value = value(i,:);

    figure(10),semilogy(value(i,:),st(i,:)/sum(st(i,:)),'color',colors(i,:)), hold on
end

figure(3),plot(Energy),hold on, plot(Energy,'+','MarkerSize',15)
save(savepath, 'files','Energy','curvature');

    
