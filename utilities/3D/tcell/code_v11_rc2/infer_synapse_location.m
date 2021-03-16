function [SaveName] = infer_synapse_location(annotation_filename, FilePath1, inferred_annotation_directory, param)
% 
% Author: Xiongtao Ruan and Long Liu
% 
% 03/06/2019 fix bug for the situation when there is no Overlay image 
%

if nargin < 4
    param.verbose = false;
    param.intermediate_location = './temp/';
    param.save_intermediate_results = false;
    param.input_mode = 'new';
end

param = process_options_structure(struct('filter_touching_tcell', false, ...
                                         'filter_touching_tcell_dist_thrsh', 10 ...
                                         ), param);

intermediate_location = param.intermediate_location;
mkdir(intermediate_location);

if ~exist(inferred_annotation_directory, 'dir')
    mkdir(inferred_annotation_directory);
end
[~, filename] = fileparts(annotation_filename);
Result_filename = ['raw_', filename, '.csv'];

% save intermediate frames of the movie
images_intermediate_path = [intermediate_location, filesep, filename, filesep];
mkdir(images_intermediate_path);
CellProcessDir = [images_intermediate_path, filesep, 'Processed/'];
mkdir(CellProcessDir);

switch param.input_mode
    case 'old'
        excel_file = csvread(annotation_filename,1,0);%csv file
        % All_StartPoint=excel_file;
        All_StartPoint_Inds = find(excel_file(:, 9) == 0 & excel_file(:, 1) ~= 0);% find all start point        
        All_StartPoint = excel_file( All_StartPoint_Inds, : );    
    case 'new'
        excel_file = csvread(annotation_filename,1,0);%csv file
        All_StartPoint=excel_file;       
end
        
GFP=GfpDeal( FilePath1);         
DIC=DICDeal( FilePath1  );
try
    Overlay=OverlayDeal(StrOverlay);
catch
    Overlay = generate_overlay_from_gfp_dic(GFP, DIC);
end
%  load GFP.mat 
%  load DIC.mat
%  load Overlay.mat
Superposition=GFP;%46 frames of GFP are processed to obtain a superimposed image of 21 layers per frame

% add code to use global segmentation for all t cell to find whether there
% is any touching tcells 
if param.filter_touching_tcell
    % segment all T cells and APC cells.
    Tcell_segment_cell = cell(numel(GFP), 1);
    Tcell_segment_loc_cell = cell(numel(GFP), 1);
    APC_segment_info_cell = cell(numel(GFP), 1);

    for i = 1 : numel(GFP)
%         GFP_i = GFP{i};
%         [loc_1, I_binary] = Tcell_segmentation(GFP_i);
%         %  figure(FrameNum);clf;imshow(image2);
%         loc = mat2cell(loc_1, ones(size(loc_1, 1), 1));
%         Tcell_segment_cell{i} = I_binary;
%         Tcell_segment_loc_cell{i} = loc;

    %     % detection of large cells
    %     % [centers, radii]=DrawBigCell(image2,15,35);
    %     DIC_i = DIC{i};
    %     [centers, radii] = APC_segmentation(DIC_i, I_binary);
    %     APC_segment_info_cell{i} = {centers, radii};
    end
end

SaveName = strcat(inferred_annotation_directory,Result_filename);%Excel path and name to save the final result
All_Need=[];

for i=1:size(All_StartPoint, 1)
% for i=9
    Processdir=strcat(CellProcessDir,num2str(i));
    if ~exist(Processdir, 'dir')
        mkdir(Processdir);
    end
	switch param.input_mode
		case 'old'
			Row_excel=All_StartPoint_Inds(i);% Row of the excel table is Row_excel-1
            Excel=[excel_file(Row_excel,2:7),excel_file(Row_excel,12)];
			Start=Excel(end);
            OnsetTime =excel_file(Row_excel, end);
            % LeftX=Excel(1);LeftY=Excel(2);RightX=Excel(3);RightY=Excel(4);
			[LeftX,LeftY,RightX,RightY]=ExchangeExcelToOur(Excel);
		case 'new'
			Row_excel=All_StartPoint(i,:);% Row of the excel table is Row_excel-1
			Start=Row_excel(end);
			OnsetTime =Start;
			LeftX=Row_excel(1);LeftY=Row_excel(2);RightX=Row_excel(3);RightY=Row_excel(4);
	end	
    Excel_lp=[LeftX,LeftY];
    Excel_rp=[RightX,RightY];
    JiaoDian_Point{Start,1}=[Excel_lp;Excel_rp];
    r=20;%Empirical radius of great circle
    Excel_cp=[(Excel_lp(1)+Excel_rp(1))/2,(Excel_lp(2)+Excel_rp(2))/2];
    width1=34;width2=60;cnt=1;
    % Small cell segmentation initialization
    [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img]=InitialSmallCell2(Superposition{Start,1},width1,Excel_cp,Processdir,Start, param);
    target_lunkuo=lunkuo_of_ImgInput1;
    lunkuo_of_ImgInput{Start,1}=lunkuo_of_ImgInput1;
    ZhixinX1{Start,1}=ZhixinX;
    ZhixinY1{Start,1}=ZhixinY;
    % Large cell segmentation initialization
    [Circle_LunKuo,Center,CircleR,Circle_Center_X,Circle_Center_Y]=InitialBigCell2(Excel_lp,Excel_rp,r,ZhixinX,ZhixinY);
    % Determine the first frame of small cells in the location of large cells
    % and distinguish the left and right point
    if isempty(JiaoDian_Point{Start,1})
        Need=[0,0,0,0,Start];
    else    
        Need=[JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),Start];%��ʼ֡excel���?    end
    end
    if param.verbose        
        Processdir_big=strcat(Processdir,'/Big/',num2str(Start));
        if ~exist(Processdir_big, 'dir')
            mkdir(Processdir_big);
        end
%         image1 = imread([file_path3,strcat(int2str(Start),'.png')],'png');
        image1=DIC{Start,1};
        figure(1);imshow(image1);hold on
        plot(target_lunkuo(:,1),target_lunkuo(:,2),'r.');
        hold on
        plot(Circle_LunKuo(:,1),Circle_LunKuo(:,2),'g.');
        hold on
        plot(JiaoDian_Point{Start,1}(1,1),JiaoDian_Point{Start,1}(1,2),'r*',JiaoDian_Point{Start,1}(2,1),JiaoDian_Point{Start,1}(2,2),'y*');
        hold on
        plot(ZhixinX,ZhixinY,'yo');
        saveas(gcf,[strcat(Processdir_big,'/'),'Result.png'],'png');
     end
     
%   [ Start, Forward_Need,Need] = No_Touch( Start,Overlay,DIC,GFP,Excel_lp, Excel_rp,OnsetTime,Need,Processdir,param);
  
    for j = Start+1:length(Superposition)
        %Small cell segmentation
        [ZhixinX,ZhixinY,lunkuo_of_ImgInput1,Img, Relative, Position]=SmallCellTracking2(Superposition{j,1},ZhixinX,ZhixinY,width1,Img,Processdir,j, param);
        if Relative
            break;
        end
        if abs(ZhixinX-ZhixinX1{j-1,1})>15 && abs(ZhixinY-ZhixinY1{j-1,1})>15
            lunkuo_of_ImgInput1=lunkuo_of_ImgInput{j-1,1};
            ZhixinX=ZhixinX1{j-1,1};
            ZhixinY=ZhixinY1{j-1,1};
        end
        lunkuo_of_ImgInput{j,1}=lunkuo_of_ImgInput1;
        ZhixinX1{j,1}=ZhixinX;
        ZhixinY1{j,1}=ZhixinY;
        %Great cell segmentation
        [Center,CircleR,Circle_LunKuo,Circle_Center_X,Circle_Center_Y]=BigCellTracking2(DIC{j,1},j,Center,CircleR,Circle_Center_X,Circle_Center_Y,width2,Superposition{j,1},Processdir, param);  
        
%         image1 = imread([file_path3,strcat(int2str(j),'.png')],'png');
        image1=DIC{j,1};
        if param.verbose
            figure(1);clf,imshow(image1);hold on
            plot(Circle_Center_X,Circle_Center_Y,'ro');
            hold on
        end
        for i_q=j-1:-1:Start
            if ~isempty(JiaoDian_Point{i_q,1})
                JiaoDian=JiaoDian_Point{i_q,1};
                break;
            end
        end
        if isempty(lunkuo_of_ImgInput1)
            for i_q=j-1:-1:Start
                if ~isempty(lunkuo_of_ImgInput{i_q,1})
                    lunkuo_of_ImgInput1=lunkuo_of_ImgInput{i_q,1};
                    ZhixinX=ZhixinX1{i_q,1};
                    ZhixinY=ZhixinY1{i_q,1};
                    break;
                end
            end
        end
        %Calculate intersection point    
        [JiaoDian_Point{j,1},cnt]=DrawPointsInLoop4(j,lunkuo_of_ImgInput1,Circle_LunKuo,Center,CircleR,cnt,JiaoDian, Superposition{j, 1}, param);
        % use similar method for infer -1/-2 to adjust the predicted
        % synapse location, so that they are not close to the centroid. 
        if ~false && ~isempty(JiaoDian_Point{j,1})
            GFP_raw = GFP{j,1};
            GFP_raw_crop = GFP_raw(round(Position(2)) + (1:size(Img, 1))-1, round(Position(1)) + (1:size(Img, 2))-1);
            I_1 = imgaussfilt(GFP_raw_crop, 1);
            thresh = graythresh(I_1);
            GFP_seg = I_1 > (thresh * 255);
            GFP_seg = imfill(GFP_seg, 'hole');
            CC  = bwconncomp(GFP_seg);
            if CC.NumObjects > 1
                % chose the largest one
                [~, max_ind] = max(cellfun(@numel, CC.PixelIdxList));
                max_PixelIdxList = CC.PixelIdxList{max_ind};
                GFP_seg = zeros(size(GFP_seg));
                GFP_seg(max_PixelIdxList) = true;
            end
            GFP_raw_seg = double(GFP_raw_crop) .* double(GFP_seg);            
            stats = regionprops(GFP_seg);
            centroid = stats.Centroid;
            JiaoDian_raw = reshape(JiaoDian_Point{j, 1}', 1, []);
            left_right_end_points = JiaoDian_raw - [Position(1:2), Position(1:2)];
            centroid_apc = [Circle_Center_X, Circle_Center_Y] - Position(1:2);
            Point = AdjustSynapse( GFP_raw_seg, GFP_seg, centroid, centroid_apc, left_right_end_points);            
            if ~isempty(Point)
                Point = Point + [Position(1:2), Position(1:2)];
                JiaoDian_Point{j,1} = [Point(1:2); Point(3:4)];
            end
        end
        
        %Determine the right and left intersection points
        if ~isempty(JiaoDian_Point{j,1})
            Need=[Need;JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),j];
            if param.verbose
                plot(JiaoDian_Point{j,1}(1,1),JiaoDian_Point{j,1}(1,2),'r*',JiaoDian_Point{j,1}(2,1),JiaoDian_Point{j,1}(2,2),'y*');
    %         saveas(gcf,[directory,int2str(j)],'png');
            end
        else
            Need=[Need;0,0,0,0,j];
        end 
        if param.verbose
            Processdir_result=strcat(Processdir,'/result');
            if ~exist(Processdir_result, 'dir')
                mkdir(Processdir_result);
            end
             saveas(gcf,[strcat(Processdir_result,'/'), int2str(j), '.png'],'png');
        end
    end
    JiaoDian_Point=[];
    Need = Interpolation_3(Need);
    Need = [Need, OnsetTime * ones(size(Need, 1), 1)];
    %确定未接触的细胞
    if size(Need,1) > 2
             Forward_Need = [];
            [ Forward_Need,Need] = No_Touch( Start,Overlay,DIC,GFP,Excel_lp, Excel_rp,OnsetTime,Need,Processdir,param);
            Forward_Need = flipud(Forward_Need);
    end
    if size(Need, 1) < 2 && strcmp(param.input_mode, 'old')
        Forward_Need = [];
    elseif size(Need, 1) < 2 && strcmp(param.input_mode, 'new')
        Forward_Need = [];
        Need = [];
    end
    
    cur_synapses_mat = [Forward_Need;Need];
	if ~isempty(cur_synapses_mat)
		cur_synapses_mat(any(cur_synapses_mat(:, 1 : 4) == 0, 2), :) = [];
		% detect whether there is touching t cells 
		if param.filter_touching_tcell
			% [is_touching] = touching_tcell_filtering_function(Tcell_segment_cell, Tcell_segment_loc_cell, [Forward_Need;Need], param);
			[is_shifting, bad_timepoint_inds] = synapse_shifting_detection_function(cur_synapses_mat, param);
			if is_shifting
				cur_synapses_mat = [];     
			else
				cur_synapses_mat(bad_timepoint_inds, :) = [];
			end
		end
	end

    All_Need=[All_Need; cur_synapses_mat; -1,-1,-1,-1,-1,-1;-1,-1,-1,-1,-1,-1];
    Need=[];
end
%Save result
% FormHead={'Lx','Ly','Rx','Ry','frame', 'onset_time'};
FormHead={'Lx','Ly','Rx','Ry','frame', 'onset_time'};
All_Need=num2cell(All_Need);
Result = cell2table(All_Need,'VariableNames',FormHead);
writetable(Result,SaveName,'WriteRowNames',true);

if ~param.save_intermediate_results
    % delete intermediate images
    rmdir(images_intermediate_path, 's');
end
% All_Need=[FormHead;All_Need];
% xlswrite(SaveName,All_Need);
end


