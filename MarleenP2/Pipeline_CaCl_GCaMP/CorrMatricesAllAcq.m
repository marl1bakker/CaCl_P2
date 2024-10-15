function CorrMatricesAllAcq(dataname, GSR)
%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

if ~exist('GSR', 'var')
    GSR = 0;
end
if GSR == 1
    addgsr = '_GSR';
else
    addgsr = '';
end

% allrois =  {'VisualROI_R', 'SensoryROI_R', 'MotorROI_R',...
%     'RetrosplenialROI_R', 'VisualROI_L','SensoryROI_L', ...
%     'MotorROI_L','RetrosplenialROI_L' };
allrois = {'VR', 'SR', 'MR', 'RR', 'VL', 'SL', 'ML', 'RL'}; % this is the order of saving in timecourses

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');

Acquisitions = {'A1', 'A2', 'A3'};
SaveFolderFig = '/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/';

if matches(dataname, 'hemoCorr_fluo.dat')
    maintitle = 'GCaMP';
else 
    maintitle = dataname(1:end-4);
end   

if ~exist([SaveFolderFig filesep maintitle], 'dir')
    mkdir([SaveFolderFig filesep maintitle])
end

%% Choose grouping
Grouping = {'Group'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);
Grouping = char(Grouping);

%% Save Av. Correlation matrixes
% Single group Corr matrix
% Start going per Acq, per group, per mouse

for indAcq = 1:size(Acquisitions, 2) % go per acuisition
    Acquisition = Acquisitions{indAcq};
    
    for index = 1:size(groups,1) % Go per group
        group = groups{index};
%         disp(group);
        
        % if exist([SaveFolderFig filesep maintitle filesep group '_' Acquisition '_corrmatrices' addgsr '.mat'], 'file')
        %     load([SaveFolderFig filesep maintitle filesep group '_' Acquisition '_corrmatrices' addgsr '.mat'], 'allcorrmatrices')
        %     eval(['AvCorrMat' group Acquisition ' = mean(allcorrmatrices, 1, ''omitnan'');'])
        %     continue
        % end

        % eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
        idx = RecordingOverview.(Grouping) == group;
        Mousegroup = RecordingOverview(idx,:);
        eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
        
        allcorrmatrices = NaN(size(Mousegroup,1), size(allrois,2), size(allrois,2));   %mice   x   roixroi (corrmatrix)
        
        for ind = 1:size(Mousegroup, 1) %go per mouse
            Mouse = Mousegroup.Mouse{ind};
            % eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
            DataFolder = [Mousegroup.(Acquisition){ind} filesep];
            SaveFolder = [Mousegroup.SaveDirectory{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
            
            if matches(Mouse, 'M23') && matches(dataname, 'hemoCorr_fluo.dat')
                continue %M23 has not enough fluorescence, so skip for fluo.
            elseif matches(Mouse, 'M14') && matches(Acquisition, 'A3') && ~matches(dataname, 'hemoCorr_fluo.dat')
                continue %M14 has really weird outliers for hbo/hbr data
            elseif matches(Mouse, 'M32') && matches(Acquisition, 'A2') && ~matches(dataname, 'hemoCorr_fluo.dat')
                continue %M32 has damaged window, hbo/hbr data looks very bizarre
            elseif matches(Mouse, 'M32') && matches(Acquisition, 'A3') && ~matches(dataname, 'hemoCorr_fluo.dat')
                continue %M32 has damaged window, hbo/hbr data looks very bizarre
            end
            
            if exist([SaveFolder 'timecourses_' dataname(1:end-4) '_centroids' addgsr '.mat'], 'file')
                load([SaveFolder 'timecourses_' dataname(1:end-4) '_centroids' addgsr '.mat'], 'AllRois');
                            
                Timecourses = cell2mat(AllRois(:,2));
                load([SaveFolder 'MovMask.mat'], 'MovMask');
                Timecourses = Timecourses .* MovMask; %to remove movement. I checked, the difference between corr with and without movement is very small, mostly 0.01 differences. checked for one mouse (M15-A1)
                load([SaveFolder 'OutlierMask.mat'], 'OutlierFrames');
                Timecourses = Timecourses .* OutlierFrames; 
                Timecourses(Timecourses == 0) = nan;
                Timecourses = Timecourses(:,1:9000);
                Timecourses(3,:) = []; % remove auditory cortex
                Timecourses(7,:) = [];

                eval([group 'N = ' group 'N+1;'])
            else
                Timecourses = NaN(8, 9000);
            end
            
            CorrMat = corr(Timecourses', 'rows','pairwise');
            
            allcorrmatrices(ind,:,:) = CorrMat;
            roi_order = AllRois(:,3);                 
            roi_order(3) = [];
            roi_order(7) = [];
            clear Timecourses CorrMat AllRois MovMask
            
        end
        
        allcorrmatrices = reshape(allcorrmatrices,[],64);
        % keep mice with nan (this is needed for stats comparing a1-a2-a3)
        % allcorrmatrices(isnan(allcorrmatrices(:,1)),:) = []; % get rid of nan mouse
                
        save([SaveFolderFig filesep maintitle filesep group '_' Acquisition '_corrmatrices' addgsr], 'allcorrmatrices', 'roi_order')
        
        AvCorrMat = mean(allcorrmatrices, 1, 'omitnan');
        eval(['AvCorrMat' group Acquisition ' = AvCorrMat;'])

    end
end

clear AvCorrMat DataFolder idx ind indAcq index Mouse Mousegroup ROItype 

%% Plot in tiledlayout
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout(3,3, 'TileSpacing', 'tight', 'Padding', 'tight');

if GSR == 0
    axlims = [0 1];
    axlimsdiff = [-0.5 0.5];
else
    axlims = [-1 1];
    axlimsdiff = [-1 1];
end

%% A1
PlotTile(AvCorrMatCaClA1, allrois,axlims)
title('CaCl', 'FontSize', 11)
ylabel('A1', 'FontWeight', 'bold', 'FontSize', 11)
PlotTile(AvCorrMatShamA1, allrois,axlims)
title('Sham', 'FontSize', 11)

Diff = AvCorrMatCaClA1 - AvCorrMatShamA1;
PlotTile(Diff, allrois, axlimsdiff);
title('Difference', 'FontSize', 11)

%% A2
PlotTile(AvCorrMatCaClA2, allrois, axlims)
ylabel('A2', 'FontWeight', 'bold', 'FontSize', 11)
PlotTile(AvCorrMatShamA2, allrois, axlims)

Diff = AvCorrMatCaClA2 - AvCorrMatShamA2;
PlotTile(Diff, allrois, axlimsdiff);

%% A3
PlotTile(AvCorrMatCaClA3, allrois, axlims)
ylabel('A3', 'FontWeight', 'bold', 'FontSize', 11)
colorbar('southoutside')
PlotTile(AvCorrMatShamA3, allrois, axlims)
colorbar('southoutside')

Diff = AvCorrMatCaClA3 - AvCorrMatShamA3;
PlotTile(Diff, allrois, axlimsdiff);
colorbar('southoutside')
 
if GSR == 0
    title(t, maintitle)
else
    title(t, [maintitle ' with GSR'])
end
f.Position = [100 70 800 800];
%% Save
% saveas(gcf, [SaveFolderFig maintitle '.tiff'], 'tiff');
% saveas(gcf, [SaveFolderFig maintitle '.eps'], 'epsc');
saveas(gcf, [SaveFolderFig maintitle addgsr '.svg'], 'svg');   
saveas(gcf, [SaveFolderFig maintitle addgsr '.png'], 'png');   
close(f)
end


function PlotTile(data, allrois, axlims)
    nexttile
    ax = gca;
    [allrois, sortindex] = sort(allrois); % alfabetisch - get way to sort
    data = reshape(data, size(allrois,2), size(allrois,2)); % alfabetisch
    data = data(sortindex, sortindex); % alfabetisch
    data = tril(data); %alfabethsich
    % data = tril(reshape(data, size(allrois,2), size(allrois,2))); %non alfabetical
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN; 

    imagesc(ax, data, 'AlphaData', ~isnan(data), axlims)
    ax.Box = 'off';
    axis image
    axis equal;
        
    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca, 'YTickLabel', ay);
%     set(gca,'YTickLabel',ay,'FontSize',15, 'FontWeight', 'bold', 'Linewidth', 2);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)
    load('/home/mbakker/P2_scripts/MarleenP2/NL.mat');
    colormap(NL)
%     colorbar
end







