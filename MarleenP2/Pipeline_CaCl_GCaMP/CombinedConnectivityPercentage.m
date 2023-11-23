function CombinedConnectivityPercentage(dataname, Acquisition, Grouping, threshold, Overwrite)

%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

if ~exist('threshold', 'var')
    threshold = 0.6;
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end

%% check if you can combine groups
% if you did the combi, you can combine it into ' group' and dont have to
% calculate new ones
% if matches(Grouping, 'Group') && ... %         ~exist([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep 'CaCl_' Acquisition '.mat'], 'file') && ...
%         exist([SaveDir '/ConnectivityPercentage/AllConnValues_CaClFemale_' dataname(1:end-4) '_' Acquisition '.mat'], 'file') && ...
%         exist([SaveDir '/ConnectivityPercentage/AllConnValues_CaClMale_' dataname(1:end-4) '_' Acquisition '.mat'], 'file')
%     
%     %cacl
%     load([SaveDir '/ConnectivityPercentage/AllConnValues_CaClFemale_' dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
%     caclfemaleconn = AllConnValues;
%     caclfemaleperc = AllPercentages;
% 
%     load([SaveDir '/ConnectivityPercentage/AllConnValues_CaClMale_' dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
%     caclmaleconn = AllConnValues;
%     caclmaleperc = AllPercentages;
%     
%     AllConnValues = [caclfemaleconn; caclmaleconn];
%     AllPercentages = [caclfemaleperc, caclmaleperc];
%     save([SaveDir '/ConnectivityPercentage/AllConnValues_CaCl_' dataname(1:end-4) '_' Acquisition '.mat'], ...
%         'AllConnValues', 'AllPercentages');
% 
%     clear AllConnValues AllPercentages cacl*
%     
%     %sham
%     load([SaveDir '/ConnectivityPercentage/AllConnValues_ShamFemale_' dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
%     shamfemaleconn = AllConnValues;
%     shamfemaleperc = AllPercentages;
% 
%     load([SaveDir '/ConnectivityPercentage/AllConnValues_ShamMale_' dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
%     shammaleconn = AllConnValues;
%     shammaleperc = AllPercentages;
%     
%     AllConnValues = [shamfemaleconn, shammaleconn];
%     AllPercentages = [shamfemaleperc, shammaleperc];
%     save([SaveDir '/ConnectivityPercentage/AllConnValues_Sham_' dataname(1:end-4) '_' Acquisition '.mat'], ...
%         'AllConnValues', 'AllPercentages');
%     
%     clear AllConnValues AllPercentages sham*
% end


%% Single group Corr matrix
% Start going per group, per mouse
for index = 1:size(groups,1) % Go per group
    
    group = groups{index};
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
    Mousegroup = RecordingOverview(idx,:);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    AllConnValues = [];
    AllPercentages = [];
    
    if exist([SaveDir '/ConnectivityPercentage/AllConnValues_' group '_' dataname(1:end-4) '_' Acquisition '.mat'], 'file') && Overwrite == 0
        load([SaveDir '/ConnectivityPercentage/AllConnValues_' group '_' dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
        
        eval(['AllConnValues_' group '= AllConnValues;']);
        eval(['AllPercentages_' group '= AllPercentages;']);
        eval([group 'N = sum(~isnan(mean(AllPercentages_' group ', 2, ''omitnan'')));'])
        
    else
        if exist([SaveDir '/ConnectivityPercentage/AllConnValues_' group dataname(1:end-4) '_' Acquisition '.mat'], 'file')
            disp('Connectivity percentage already done, OVERWRITING MAT FILES')
        end
        
        for ind = 1:size(Mousegroup, 1) %go per mouse
            Mouse = Mousegroup.Mouse{ind};
            eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
            SaveFolder = [Mousegroup.SaveDirectory{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
            
            [ConnValuesMouse, nrpixels] = SingleSubjectConnectivityValues(SaveFolder, dataname, threshold);
            PercentagesMouse = ConnValuesMouse/nrpixels;
            
            AllConnValues (ind,:) = ConnValuesMouse;
            AllPercentages(ind,:) = PercentagesMouse;
        end
        
        save([SaveDir '/ConnectivityPercentage/AllConnValues_' group dataname(1:end-4) '_' Acquisition '.mat'], 'AllConnValues', 'AllPercentages');
        
        eval(['AllConnValues_' group '= AllConnValues;']);
        eval(['AllPercentages_' group '= AllPercentages;']);
        eval([group 'N = sum(~isnan(mean(AllPercentages_' group ', 2, ''omitnan'')));'])
        
    end
end

clear AllConnValues AllPercentages idx index indx Possibilities 


%% Make Table
overviewtable = table;

for indgroup = 1:size(groups,1)
    group = groups{indgroup};
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
    Mousegroup = RecordingOverview(idx,:);

    eval(['toplot = AllPercentages_' group ';']) 
    toplot = toplot';
    toplot(3,:) = []; %take out auditory
    toplot(7,:) = []; %toplot is roi x mice
    %     tablegroup = table;
    labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};
    
    for indmouse = 1:size(Mousegroup, 1)
        if matches(Mousegroup.Mouse{indmouse}, 'M23') && matches(dataname, 'hemoCorr_fluo.dat')
            eval([group 'N = ' group 'N - 1;'])
            continue
        end
        tablemouse = table;
        tablemouse.Mouse = cellstr(repmat(Mousegroup.Mouse{indmouse}, size(labels,2),1));
%         tablemouse.Acquisition = cellstr(repmat(Acquisition, size(labels,2),1));
        tablemouse.ROI = labels';
        tablemouse.Group = cellstr(repmat(group, size(labels, 2),1));
        tablemouse.PercentageCorr = toplot(:,indmouse);
        overviewtable = [overviewtable; tablemouse];
    end
end

overviewtable.Group = categorical(overviewtable.Group);
overviewtable.ROI = categorical(overviewtable.ROI);
overviewtable.PercentageCorr = overviewtable.PercentageCorr.*100;

%% start plotting
f = figure('InvertHardcopy','off','Color',[1 1 1]);

% if sexgrouping == 1 % if you're grouping by sex
%     caclind = contains(string(overviewtable.Group), 'CaCl');
%     cacl = overviewtable(caclind,:);
%     shamind = contains(string(overviewtable.Group), 'Sham');
%     sham = overviewtable(shamind,:);
%     
%     b = tiledlayout(2,1);
%     axes = nexttile(b);
%     boxchart(cacl.ROI, cacl.Fluctuations, 'GroupByColor', cacl.Group,...
%         'LineWidth', 2);
%     
%     %Make pretty
%     legend( 'Location', 'northeast', 'NumColumns', 2)
%     ylabel('Fluctuation (std)')
%     if matches(dataname, 'hemoCorr_fluo.dat')
%         ylim([0 0.02])
%     else
%         ylim([0 1])
%     end
%     set(axes,'FontSize',20,'FontWeight','bold','LineWidth',2);
% 
%     % next tile
%     axes = nexttile(b);
%     boxchart(sham.ROI, sham.Fluctuations, 'GroupByColor', sham.Group,...
%         'LineWidth', 2);
%     legend( 'Location', 'northeast', 'NumColumns', 2)
%     
%     %Make pretty
%     legend( 'Location', 'northeast', 'NumColumns', 2)
%     ylabel('Fluctuation (std)')
%     if matches(dataname, 'hemoCorr_fluo.dat')
%         ylim([0 0.02])
%     else
%         ylim([0 1])
%     end
%     set(axes,'FontSize',20,'FontWeight','bold','LineWidth',2);
%     
%     Grouping = 'Sex';
%     title(b, ['Fluctuation ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')

% else

% without scatters:
% b = boxchart(overviewtable.ROI, overviewtable.PercentageCorr, 'GroupByColor', overviewtable.Group,...
%     'LineWidth', 2);

% overlay scatter plot
overviewtable.idx = grp2idx(overviewtable.ROI); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = sort(labels);
b = boxchart(overviewtable.idx, overviewtable.PercentageCorr, 'GroupByColor', overviewtable.Group,...
    'LineWidth', 2);

for ind = 1:size(groups,1)
    b(ind).MarkerStyle = 'none'; %dont display outliers because we will do scatter that will show them
end
hold on

% overlay the scatter plots
xaxisplacement = 0.375; %start position

for indroi = 1:length(labels)
    currentROI = overviewtable.ROI == labels{indroi};
    
    for indgroup = 1:size(groups,1)
        group = groups{indgroup};
        currentgroup = overviewtable.Group == group;
        currentindex = currentgroup.*currentROI;
        overviewtable.idx(currentindex == 1) = xaxisplacement + 0.25*indgroup;
    end
    
    xaxisplacement = xaxisplacement + 1;
end
clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot
 
hs = scatter(overviewtable.idx, overviewtable.PercentageCorr, 'filled', 'jitter','on','JitterAmount',0.02);
hs.MarkerFaceColor = [0 0 0];
hs.MarkerFaceAlpha = 0.3;
xticks(1:length(labels));
xticklabels(labels)
xlim([0.2 length(labels)+0.7])

%Make pretty
axes1 = b.Parent;
hold(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
legend(b, 'Location', 'northeast', 'NumColumns', 2)
ylim([0 100])
ylabel(['% pixels with correlation > ' num2str(threshold)])
xlabel('Region of Interest')

b(1).SeriesIndex = 7;
b(2).SeriesIndex = 1;
if size(b,1)>2
    b(3).SeriesIndex = 2;
    b(4).SeriesIndex = 6;
end

title(['Connectivity ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
subtitleN = [];
for indgroup = 1:size(groups, 1)
    eval(['ngroup = num2str(' groups{indgroup} 'N);'])
    subtitleN = [subtitleN groups{indgroup} ' N = ' ngroup ' -- '];
end
subtitle(subtitleN(1:end-4), 'FontSize', 10)

f.Position = [10 10 1500 1000]; %for size of screen before saving

if ~exist([SaveDir '/ConnectivityPercentage/' dataname(1:end-4)], 'dir')
    mkdir([SaveDir '/ConnectivityPercentage/' dataname(1:end-4)])
end

saveas(gcf, [SaveDir '/ConnectivityPercentage/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/ConnectivityPercentage/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.eps'], 'epsc');

close(f)

end