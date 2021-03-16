% function [cPrecision,cRecall] = evalseg( manual, water)
function [cPrecision,cRecall,pPrecision,pRecall,counts] = evalseg( manual, water)
% clear
% load tmpo1

MM = uint16(bwlabel(logical(water),4));
HM = manual;

overlapHM = double(uint16(MM>0).*manual);
overlapMM = bwlabel(logical(MM),4).*(manual>0);

statsHM = [];
nHM = double(max(HM(:)));

propsHM = regionprops(HM,'Area');
for i=1:nHM

    statsHM(i).pixel_area = propsHM(i).Area;
    tmp = (HM==i).*overlapMM;
    
    u = unique(tmp(tmp>0));
    statsHM(i).mach_areasHM = zeros(1,length(u));
    for j=1:length(u)
        statsHM(i).mach_areasHM(j) = sum(tmp(:)==u(j));
    end
    [ma ind] = max(statsHM(i).mach_areasHM);
    statsHM(i).overlapHM = ma / statsHM(i).pixel_area;
    statsHM(i).mach_ind = u(ind);
end

mach_ind = zeros(1,nHM);
overlap_val = zeros(1,nHM);
for i=1:nHM
    if isempty(statsHM(i).mach_ind)
        mach_ind(i) = 0;
    else
        mach_ind(i) = statsHM(i).mach_ind;
    end
    if isempty(statsHM(i).overlapHM)
        overlap_val(i) = 0;
    else
        overlap_val(i) = statsHM(i).overlapHM;
    end
end

[ccc bbb] = hist(mach_ind(mach_ind>0),1:max(mach_ind));
indccc = find(ccc==1);

nMM = double(max(MM(:)));

thr = 0.75;
cTP = sum(overlap_val>thr & ismember(mach_ind,indccc));
cFP = nMM-cTP;
cFN = nHM - cTP;

cRecall = cTP/(nHM);
cPrecision = cTP/(cTP+cFP);


% return

statsMM = [];
for i=1:nMM
    statsMM(i).pixel_area = sum(MM(:)==i);
    tmp = (MM==i).*overlapHM;
    
    u = unique(tmp);
    u(u==0) = [];
    statsMM(i).hand_areasMM = zeros(1,length(u));
    for j=1:length(u)
        statsMM(i).hand_areasMM(j) = sum(tmp(:)==u(j));
    end
    [ma ind] = max(statsMM(i).hand_areasMM);
    statsMM(i).overlapMM = ma / statsMM(i).pixel_area;
    statsMM(i).hand_ind = u(ind);
    
    if statsHM(u(ind)).mach_ind==i
        statsMM(i).npixel_area = sum(statsMM(i).hand_areasMM) - ma;
    else
        statsMM(i).npixel_area = sum(statsMM(i).hand_areasMM);
    end
end

pTP = 0;
for i=1:nHM
    if isempty(statsHM(i).overlapHM)
        pTP = pTP + 0;
    else
        pTP = pTP + statsHM(i).overlapHM*statsHM(i).pixel_area;
    end
end

pFP = 0;
for i=1:nMM
    pFP = pFP + statsMM(i).npixel_area;
end

pFN = (sum(MM(:)>0)) - pTP;

pRecall = pTP/(sum(MM(:)>0));
pPrecision = pTP/(pTP+pFP);


counts = [cTP cFP cFN pTP pFP pFN];















% % % function [cPrecision,cRecall] = evalseg( manual, water)
% % function [cPrecision,cRecall,pPrecision,pRecall,counts] = evalseg( manual, water)
% % % clear
% % % load tmpo1
% % 
% % MM = uint16(bwlabel(logical(water),4));
% % HM = manual;
% % 
% % overlapHM = double(uint16(MM>0).*manual);
% % overlapMM = bwlabel(logical(MM),4).*(manual>0);
% % 
% % statsHM = [];
% % nHM = double(max(HM(:)));
% % 
% % propsHM = regionprops(HM,'Area');
% % for i=1:nHM
% % 
% %     statsHM(i).pixel_area = propsHM(i).Area;
% %     tmp = (HM==i).*overlapMM;
% %     
% %     u = unique(tmp(tmp>0));
% %     statsHM(i).mach_areasHM = zeros(1,length(u));
% %     for j=1:length(u)
% %         statsHM(i).mach_areasHM(j) = sum(tmp(:)==u(j));
% %     end
% %     [ma ind] = max(statsHM(i).mach_areasHM);
% %     statsHM(i).overlapHM = ma / statsHM(i).pixel_area;
% %     statsHM(i).mach_ind = u(ind);
% % end
% % 
% % mach_ind = zeros(1,nHM);
% % overlap_val = zeros(1,nHM);
% % for i=1:nHM
% %     if isempty(statsHM(i).mach_ind)
% %         mach_ind(i) = 0;
% %     else
% %         mach_ind(i) = statsHM(i).mach_ind;
% %     end
% %     if isempty(statsHM(i).overlapHM)
% %         overlap_val(i) = 0;
% %     else
% %         overlap_val(i) = statsHM(i).overlapHM;
% %     end
% % end
% % 
% % [ccc bbb] = hist(mach_ind(mach_ind>0),1:max(mach_ind));
% % indccc = find(ccc==1);
% % 
% % nMM = double(max(MM(:)));
% % 
% % thr = 0.75;
% % cTP = sum(overlap_val>thr & ismember(mach_ind,indccc));
% % cFP = nMM-cTP;
% % 
% % cRecall = cTP/(nHM);
% % cPrecision = cTP/(cTP+cFP);
% % 
% % cFN = nHM - cTP;
% % 
% % % return
% % 
% % statsMM = [];
% % for i=1:nMM
% %     statsMM(i).pixel_area = sum(MM(:)==i);
% %     tmp = (MM==i).*overlapHM;
% %     
% %     u = unique(tmp);
% %     u(u==0) = [];
% %     statsMM(i).hand_areasMM = zeros(1,length(u));
% %     for j=1:length(u)
% %         statsMM(i).hand_areasMM(j) = sum(tmp(:)==u(j));
% %     end
% %     [ma ind] = max(statsMM(i).hand_areasMM);
% %     statsMM(i).overlapMM = ma / statsMM(i).pixel_area;
% %     statsMM(i).hand_ind = u(ind);
% %     
% %     if statsHM(u(ind)).mach_ind==i
% %         statsMM(i).npixel_area = sum(statsMM(i).hand_areasMM) - ma;
% %     else
% %         statsMM(i).npixel_area = sum(statsMM(i).hand_areasMM);
% %     end
% % end
% % 
% % pTP = 0;
% % for i=1:nHM
% %     if isempty(statsHM(i).overlapHM)
% %         pTP = pTP + 0;
% %     else
% %         pTP = pTP + statsHM(i).overlapHM*statsHM(i).pixel_area;
% %     end
% % end
% % 
% % pFP = 0;
% % for i=1:nMM
% %     pFP = pFP + statsMM(i).npixel_area;
% % end
% % 
% % pRecall = pTP/(sum(MM(:)>0));
% % pPrecision = pTP/(pTP+pFP);
% % 
% % pFN = (sum(MM(:)>0)) - pTP;
% % 
% % counts = [cTP cFP cFN pTP pFP pFN];











