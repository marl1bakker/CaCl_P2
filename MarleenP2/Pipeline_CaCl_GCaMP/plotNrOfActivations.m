function plotNrOfActivations(Acq)

% get data
try
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/NrOfActivations.mat', 'ActsTable')
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches_random.mat', 'PercTable')
    PercTable_random = PercTable;
    load('/media/mbakker/GDrive/P2/GCaMP/NVC/PercentageMatches.mat', 'PercTable')
catch
    disp('Tables not found, run NumberOfActivations first!')
    return
end


if ~exist('Acq', 'var')
    Acq = 'A1';
end

Tables.Activities = ActsTable(matches(ActsTable.Acq, Acq),:);
Tables.Percentages = PercTable(matches(PercTable.Acq, Acq),:);
Tables.Perc_random = PercTable_random(matches(PercTable_random.Acq, Acq),:);

clear PercTable_random PercTable ActsTable

%% tables in right format
tablesorts = fieldnames(Tables);
ynames = {'Act', 'Perc', 'Perc'};

for indtables = 1:size(tablesorts,1)
    currenttable = Tables.(tablesorts{indtables});


    temp = table;
    for ind = 1:size(currenttable,1) % go per mouse
        mousetable = rows2vars(currenttable(ind,5:15)); %hardcoded
        mousetable.Properties.VariableNames(1) = "ROIs";
        mousetable.Properties.VariableNames(2) = ynames(indtables);
        mousetable = [repmat(currenttable(ind,1:4), size(mousetable,1),1), mousetable]; %hardcoded 1:4

        temp = [temp; mousetable];
    end

    % get percentages 0-100 instead of 0-1
    if matches(ynames{indtables}, 'Perc')
        temp.Perc = temp.Perc.*100;

    % elseif matches(ynames{indtables}, 'Perc_r')
    %     temp.Perc_r = temp.Perc_r.*100;
    end

    temp(contains(temp.ROIs, 'Auditory'),:) = [];
    wholebrain.(tablesorts{indtables}) = temp(matches(temp.ROIs, 'WholeBrain'),:);
    temp(matches(temp.ROIs, 'WholeBrain'),:) = [];
    % Save table
    Tables.(tablesorts{indtables}) = temp;
    clear temp mousetable
end


% get averages
% disp(tablesorts{indtables})
tablesorts = {'Percentages', 'Perc_random'};
for indpercentages = 1:length(tablesorts)
    disp(tablesorts{indpercentages})

    wbhbo = wholebrain.(tablesorts{indpercentages})(matches(wholebrain.(tablesorts{indpercentages}).DetectedOn, 'HbO'),:);
    disp(['hbo: ' num2str(mean(wbhbo.Perc))])
    wbfluo = wholebrain.(tablesorts{indpercentages})(matches(wholebrain.(tablesorts{indpercentages}).DetectedOn, 'hemoCorr_fluo'),:);
    disp(['fluo: ' num2str(mean(wbfluo.Perc))])
    clear wbfluo wbhbo

    shamhbo = wholebrain.(tablesorts{indpercentages})(wholebrain.(tablesorts{indpercentages}).Group == 'Sham',:);
    shamhbo = shamhbo(matches(shamhbo.DetectedOn, 'HbO'),:);
    disp(['sham hbo: ' num2str(mean(shamhbo.Perc))])
    clear shamhbo

    caclhbo = wholebrain.(tablesorts{indpercentages})(wholebrain.(tablesorts{indpercentages}).Group == 'CaCl',:);
    caclhbo = caclhbo(matches(caclhbo.DetectedOn, 'HbO'),:);
    disp(['cacl hbo: ' num2str(mean(caclhbo.Perc))])
    clear caclhbo

    shamfluo = wholebrain.(tablesorts{indpercentages})(wholebrain.(tablesorts{indpercentages}).Group == 'Sham',:);
    shamfluo = shamfluo(matches(shamfluo.DetectedOn, 'hemoCorr_fluo'),:);
    disp(['sham fluo: ' num2str(mean(shamfluo.Perc))])
    clear shamfluo

    caclfluo = wholebrain.(tablesorts{indpercentages})(wholebrain.(tablesorts{indpercentages}).Group == 'CaCl',:);
    caclfluo = caclfluo(matches(caclfluo.DetectedOn, 'hemoCorr_fluo'),:);
    disp(['cacl fluo: ' num2str(mean(caclfluo.Perc))])
    clear caclfluo

end


%% start plotting - Per roi
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout('flow');

% Activations
nexttile
temp_table = Tables.Activities(matches(Tables.Activities.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Act', 'Detected Activations',[0 2200000]);
title('GCaMP-detected')
subtitle('Number of activations')

nexttile
temp_table = Tables.Activities(matches(Tables.Activities.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Act', 'Detected Activations',[0 2200000]);
title('HbO-detected')
subtitle('Number of activations')

% Percentages
nexttile
temp_table = Tables.Percentages(matches(Tables.Percentages.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('GCaMP-detected')
subtitle(["Detected GCaMP activations followed", "by detected HbO activations"])

nexttile
temp_table = Tables.Percentages(matches(Tables.Percentages.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('HbO-detected')
subtitle(["Detected HbO activations preceded", "by detected GCaMP activations"])

% Random percentages
nexttile
temp_table = Tables.Perc_random(matches(Tables.Perc_random.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('Random - GCaMP')
subtitle(["Randomly chosen GCaMP points followed", "by detected HbO activations"])

nexttile
temp_table = Tables.Perc_random(matches(Tables.Perc_random.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('Random - HbO')
subtitle(["Randomly chosen HbO points preceded", "by detected GCaMP activations"])

leg = legend({'CaCl', 'Sham'}, 'Orientation', 'Horizontal');
leg.Layout.Tile = 'south';

f.Position = [10 10 1200 1000];
pause(1)
saveas(f, ['/media/mbakker/GDrive/P2/GCaMP/NVC/NrofActs/Nr_and_matches_Activations_' Acq '.svg'], 'svg')

close(f)

%% plot whole brain
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout('flow');

% Activations
nexttile
temp_table = wholebrain.Activities(matches(wholebrain.Activities.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Act', 'Detected Activations', [0 15000000]);
title('GCaMP-detected')
subtitle('Number of activations')

nexttile
temp_table = wholebrain.Activities(matches(wholebrain.Activities.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Act', 'Detected Activations', [0 15000000]);
title('HbO-detected')
subtitle('Number of activations')

% Percentages
nexttile
temp_table = wholebrain.Percentages(matches(wholebrain.Percentages.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('GCaMP-detected')
subtitle(["Detected GCaMP activations followed", "by detected HbO activations"])

nexttile
temp_table = wholebrain.Percentages(matches(wholebrain.Percentages.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('HbO-detected')
subtitle(["Detected HbO activations preceded", "by detected GCaMP activations"])

% Random Percentages
nexttile
temp_table = wholebrain.Perc_random(matches(wholebrain.Perc_random.DetectedOn, 'hemoCorr_fluo'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('Random - GCaMP')
subtitle(["Randomly chosen GCaMP points followed", "by detected HbO activations"])

nexttile
temp_table = wholebrain.Perc_random(matches(wholebrain.Perc_random.DetectedOn, 'HbO'),:);
MakeBoxplotTile(temp_table, 'ROIs', 'Perc', 'Percentage', [0 100]);
title('Random - HbO')
subtitle(["Randomly chosen HbO points preceded", "by detected GCaMP activations"])

f.Position = [10 10 1200 1000];
pause(1)
saveas(f, ['/media/mbakker/GDrive/P2/GCaMP/NVC/NrofActs/Nr_and_matches_Activations_WholeBrain_' Acq '.svg'], 'svg')

close(f)


end



function MakeBoxplotTile(overviewtable, x, y, y_label, ylimvalues)

if ~exist('x','var')
    x = 'ROI';
end

%name the columns in your table that you want to use x and y
temp = find(matches(overviewtable.Properties.VariableNames, x));
overviewtable.Properties.VariableNames{temp} = 'x';
overviewtable.x = categorical(overviewtable.x);
overviewtable.idx = grp2idx(overviewtable.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(overviewtable.x))';

if ~exist('y', 'var') && sum(contains(overviewtable.Properties.VariableNames, 'y')) == 0
    disp('Cannot make boxplot tile without Y variable')
    return
end

temp = find(matches(overviewtable.Properties.VariableNames, y));
overviewtable.Properties.VariableNames{temp} = 'y';


% BOXPLOT
%dont display outliers because we will do scatter that will show them
b = boxchart(overviewtable.idx, overviewtable.y, 'GroupByColor', overviewtable.Group, 'LineWidth', 2, 'MarkerStyle', 'none');
hold on

% SCATTER
groups = cellstr(unique(overviewtable.Group));
Grouping = 'Group';
xaxisstep = 1/size(groups,1);
xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;

for indroi = 1:length(labels)
    currentROI = overviewtable.x == labels{indroi}; %is called currentroi because it's usually grouped by ROI, but can be something else

    for indgroup = 1:size(groups,1)
        currentgroup = overviewtable.(Grouping) == groups{indgroup};
        currentindex = currentgroup.*currentROI;
        overviewtable.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
    end

    xaxisplacement = xaxisplacement + 1;
end
clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot

% Scatter per group so you can vary colours
if matches(Grouping, 'Group') % hardcoded
    cacl = overviewtable(overviewtable.(Grouping)=='CaCl',:);
    hs = scatter(cacl.idx, cacl.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0.6350 0.0780 0.1840];
    hs.MarkerFaceAlpha = 0.5;
    clear hs cacl

    sham = overviewtable(overviewtable.(Grouping)=='Sham',:);
    hs = scatter(sham.idx, sham.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0 0.4470 0.7410];
    hs.MarkerFaceAlpha = 0.5;
    clear hs sham
else
    hs = scatter(overviewtable.idx, overviewtable.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.5;
end

xticks(1:length(labels));
% xticklabels(labels);
if length(labels)>1
    xticklabels({'ML','MR','RL','RR','SL','SR','VL','VR'}) %HARDCODED
else
    xticklabels({'Whole Brain'})
end
xlim([0.2 length(labels)+0.7])

axes1 = b.Parent;
hold(axes1,'on');
set(axes1, 'FontSize', 15, 'LineWidth', 2)
% set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
% legend(b, 'Location', 'northeast', 'NumColumns', 2);
if exist('ylimvalues', 'var')
    ylim(ylimvalues);
end

b(1).SeriesIndex = 7;
b(2).SeriesIndex = 1;
if size(b,1)>2
    b(3).SeriesIndex = 2;
    b(4).SeriesIndex = 6;
end

ylabel(y_label)
end
