function CombinedGCaMPFluctuations_dotplot(Acquisition, dataname)

% set up
load('/home/mbakker/P2_scripts/MarleenP2/LargerAreaAtlas.mat', 'Map')
Map(isnan(Map)) = 11;
SaveDir = '/media/mbakker/GDrive/P2/GCaMP/';

if ~exist([SaveDir 'Fluctuations/FluctTable_' dataname '.mat'], 'file')
    CombinedGCaMPFluctuations(Acquisition, dataname)
end
load([SaveDir 'Fluctuations/FluctTable_' dataname '.mat'], 'overviewtable');
overviewtable = overviewtable(overviewtable.Acquisition == Acquisition,:);

%% get good table
PlottingTable = table();

% Get right locations on the map
DotLocationsLeft = [[150,700]; [200,700]; [250,700]; [300,700];... % Vis L
    [150,460]; [200,460]; [250,460]; [300,460];... % Sen L
    [0,0];[0,0];[0,0];[0,0];... % Aud L
    [250,220]; [300,220]; [350,220]; [400,220];... % Mot L
    [430,500]; [390,560]; [430,600]; [390, 670]]; 
DotLocationsRight = [464+(464-DotLocationsLeft(:,1)), ...
    DotLocationsLeft(:,2)]; %464 is halfway point
DotLocationsRight(9:12,1) = 0; %aud R
PlottingTable.DotLocations = [DotLocationsLeft; DotLocationsRight];

% dotlocations in order: 
PlottingTable.ROI = [repmat('Vis-L',4,1); repmat('Sen-L',4,1); repmat('Aud-L',4,1);...
    repmat('Mot-L',4,1); repmat('Ret-L',4,1);...
    repmat('Vis-R',4,1); repmat('Sen-R',4,1); repmat('Aud-R',4,1);...
    repmat('Mot-R',4,1); repmat('Ret-R',4,1)];
    
Groups = {'CaClFemale';'CaClMale';'ShamFemale';'ShamMale'};
PlottingTable.Groups = repmat(Groups,10,1);
PlottingTable.Std = NaN(40,1);

for index = 1:size(PlottingTable, 1)
    CurrentROI = PlottingTable.ROI(index,:);
    CurrentGroup = PlottingTable.Groups{index,:};
    
    ROIgroupav = overviewtable(overviewtable.Combi == CurrentGroup,:);
    ROIgroupav = ROIgroupav(ROIgroupav.ROI == CurrentROI,:);
    
    if size(ROIgroupav,1) < 1
        continue % no group with that name (auditory)
    end

    PlottingTable.Std(index) = mean(ROIgroupav.StdAct, 'all','omitnan');
end

clear index CurrentROI CurrentGroup ROIgroupav DotLocations*
PlottingTable.Groups = categorical(PlottingTable.Groups);

%% plot
f = figure('InvertHardcopy','off','Color',[1 1 1]);
imcontour(Map, 'Color', 'k')
hold on
colours = [[0.6350 0.0780 0.1840];... %caclfemale
    [0 0.4470 0.7410];... %caclmale
    [0.8500 0.3250 0.0980];... %shamfemale
    [0.3010 0.7450 0.9330]]; %shammale

for ind = 1:size(Groups)
    CurrentGroup = PlottingTable(PlottingTable.Groups == Groups{ind},:);
    
    for ind2 = 1:size(CurrentGroup,1)
        if CurrentGroup.DotLocations(ind2,1) == 0
            continue
        end
        
        p = nsidedpoly(1000, 'Center', CurrentGroup.DotLocations(ind2,:), ...
            'Radius', CurrentGroup.Std(ind2,:)*2000);
        plot(p, 'FaceColor', colours(ind,:), 'FaceAlpha', 1);
    end

end

    
%only does lines, not filled circle:
% for ind = 1:size(Groups)
%     CurrentGroup = PlottingTable(PlottingTable.Groups == Groups{ind},:);
%     viscircles(CurrentGroup.DotLocations, CurrentGroup.Std*1000, 'Color', colours(ind,:));
% end

title(['Fluctuations ' Acquisition]);
subtitle('Radius = std*2000');
saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/Fluctuations/DotPlot/' Acquisition '.tiff'], 'tiff');
saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/Fluctuations/DotPlot/' Acquisition '.eps'], 'epsc');
close(f)

end

    
    
    
    