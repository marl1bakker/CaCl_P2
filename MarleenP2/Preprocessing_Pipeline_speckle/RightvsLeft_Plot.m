function RightvsLeft_Plot(datatype)

if ~exist('datatype', 'var')
    datatype = 'speckle';
end

Acquisitions = {'A1', 'A2', 'A3'};

%% get table
if matches(datatype, 'speckle')
    ratiotablepath = '/media/mbakker/GDrive/P2/Speckle/LR_ratio/LR_ratio_table.mat';
    SaveDir = '/media/mbakker/GDrive/P2/Speckle/LR_ratio/';
    ylimvalues = [0.8 1.3];
else
    ratiotablepath = ['/media/mbakker/GDrive/P2/GCaMP/LR_ratio/LR_ratio_table_' datatype '.mat'];
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP/LR_ratio/';
    ylimvalues = [0.99 1.01];
end

load(ratiotablepath, 'LR_ratio_table');
LR_ratio_table(matches(LR_ratio_table.ROI, 'Auditory'),:) = [];
LR_ratio_table = LR_ratio_table(~ismissing(LR_ratio_table.Mouse),:);
LR_ratio_table.ROI = categorical(LR_ratio_table.ROI);

% sanity check:
% boxchart(LR_ratio_table.ROI, LR_ratio_table.Ratio, 'GroupByColor', LR_ratio_table.Group)

Grouping = {'Group'};
% [groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);
groups = {'CaCl', 'Sham'}';

%% go per acquisition
%% Plot in tiledlayout
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout(1,3, 'TileSpacing', 'tight', 'Padding', 'tight');

for indacq = 1:size(Acquisitions, 2)
    TableAcq = LR_ratio_table(matches(LR_ratio_table.Acq, Acquisitions{indacq}),:);
    PlotTile(Grouping, groups, TableAcq, ylimvalues, 'ROI', 'Ratio');
end
%% A1
% PlotTile(AvCorrMatCaClA1, allrois,[0 1])
% title('CaCl', 'FontSize', 11)
% ylabel('A1', 'FontWeight', 'bold', 'FontSize', 11)



    % [f, t] = MakeBoxplot(Grouping, groups, LR_ratio_table, ylimvalues, 'ROI', 'Ratio', 1);
% f.Position = [20 20 800 1000];
% title(['Right/Left ratio ' datatype])
% legend(t, 'Location', 'northwest', 'NumColumns', 2);



% save
saveas(f, [SaveDir 'RightVsLeft_allroi_' datatype '.tiff'], 'tiff');
% saveas(f, [SaveDir 'RightVsLeft_allroi.tiff'], 'tiff');
% saveas(f, [SaveDir 'RightVsLeft_allroi.eps'], 'epsc');
close(f)

end

function PlotTile(Grouping, groups, overviewtable, ylimvalues, x, y)
nexttile

Grouping = char(Grouping);

%name the columns in your table that you want to use x and y
temp = find(matches(overviewtable.Properties.VariableNames, x));
overviewtable.Properties.VariableNames{temp} = 'x';
overviewtable.x = categorical(overviewtable.x);

if exist('y', 'var')
    temp = find(matches(overviewtable.Properties.VariableNames, y));
    overviewtable.Properties.VariableNames{temp} = 'y';
end

overviewtable.idx = grp2idx(overviewtable.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(overviewtable.x))';


% BOXPLOT
%dont display outliers because we will do scatter that will show them
eval(['b = boxchart(overviewtable.idx, overviewtable.y, ''GroupByColor'', overviewtable.' Grouping ', ''LineWidth'', 2, ''MarkerStyle'', ''none'');'])
hold on

% SCATTER
xaxisstep = 1/size(groups,1);
xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;

for indroi = 1:length(labels)
    currentROI = overviewtable.x == labels{indroi}; %is called currentroi because it's usually grouped by ROI, but can be something else

    for indgroup = 1:size(groups,1)
        eval(['currentgroup = overviewtable.' Grouping ' == groups{indgroup};']);
        currentindex = currentgroup.*currentROI;
        overviewtable.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
    end

    xaxisplacement = xaxisplacement + 1;
end
clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot

hs = scatter(overviewtable.idx, overviewtable.y, 70, 'filled', 'jitter','on','JitterAmount',0.02);

% MAKE PRETTY
hs.MarkerFaceColor = [0 0 0];
hs.MarkerFaceAlpha = 0.3;
xticks(1:length(labels));
xticklabels(labels);
xlim([0.2 length(labels)+0.7])

axes1 = b.Parent;
hold(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
legend(b, 'Location', 'northeast', 'NumColumns', 2);
%     ylim([0 100]);
ylim(ylimvalues);

b(1).SeriesIndex = 7;
b(2).SeriesIndex = 1;
if size(b,1)>2
    b(3).SeriesIndex = 2;
    b(4).SeriesIndex = 6;
end

%     title(['Connectivity ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
subtitleN = [];
for indgroup = 1:size(groups, 1)
    eval(['ngroup = sum(overviewtable.' Grouping ' == groups{indgroup});'])
    ngroup = ngroup/length(labels); %because you have a row for each ROI and you count them all. divide by 8 and you have the nr of mice
    subtitleN = [subtitleN groups{indgroup} ' N = ' num2str(ngroup) ' -- '];
end
subtitle(subtitleN(1:end-4), 'FontSize', 15)

end