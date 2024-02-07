%% make boxplot
% Grouping is usually {'Group', 'Sex'}

function BoxplotSeedCorr(Grouping, dataname, Acquisition)

%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end

clear Possibilities indx

%% Calculate correlation matrices
for indGroup = 1:size(groups,1) % Go per group
    
    group = groups{indGroup};
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
    Mousegroup = RecordingOverview(idx,:);
    Mousegroup.CorrVal = NaN(size(Mousegroup,1),1);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    
    for indMouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indMouse};
        eval(['DataFolder = [Mousegroup.' Acquisition '{indMouse} filesep];']);
        SaveFolder = [Mousegroup.SaveDirectory{indMouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
        if exist([SaveFolder 'timecourses_' dataname '.mat'], 'file')
            load([SaveFolder 'timecourses_' dataname '.mat'], 'AllRois');
            
            Timecourses = cell2mat(AllRois(:,2));
            Timecourses = Timecourses(:,1:9000);
            Timecourses(3,:) = []; %take out auditory
            Timecourses(7,:) = [];
            
            eval([group 'N = ' group 'N+1;'])
        else
            Timecourses = NaN(8, 9000);
        end
        
        CorrMat = corr(Timecourses');
        CorrMat = tril(CorrMat);
        CorrMat(CorrMat == 0 ) = NaN;
        CorrMat(CorrMat == 1 ) = NaN;
        
        Mousegroup.CorrVal(indMouse) = mean(CorrMat, 'all', 'omitnan');
                
        clear Timecourses CorrMat AllRois 
        
    end % of mice
    eval([group ' = Mousegroup;']);

end % of groups

%% go plot
% Make table T to work with
eval(['T =' groups{1} ';'])
for indgroup = 2:size(groups,1)
    eval(['T = [T;' groups{indgroup} '];'])
end

% eval(

% % Plot boxplot, grouped by sex
% f = figure('InvertHardcopy','off','Color',[1 1 1]);    
% boxchart(T.Group, T.CorrVal, 'GroupByColor', T.Sex, ...
%     'MarkerStyle', '.', 'LineWidth', 2, 'BoxWidth', 0.75);

% Plot boxplot with individual datapoints
f = figure('InvertHardcopy','off','Color',[1 1 1]);    
eval(['b = boxchart(T.' Grouping ', T.CorrVal, ''MarkerStyle'', ''.'', ''LineWidth'', 2, ''BoxWidth'', 0.75);'])
ax = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax,'FontSize', 18, 'FontWeight', 'bold', 'Linewidth', 2);
xlabel('Experimental group')
ylabel('Correlation')
ylim([0.5 1])
f.Position = [10 10 500 1000]; 

% Plot individual datapoints:
hold on
eval(['s = scatter(T.' Grouping ', T.CorrVal, ''filled'', ''Jitter'', 1);'])
s.SizeData = 50;

% title('Average correlation value, all seeds')
saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/BoxplotSeedCorrelation/' dataname '_' Grouping '_' Acquisition '.eps'], 'epsc');
saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/BoxplotSeedCorrelation/' dataname '_' Grouping '_' Acquisition '.tiff'], 'tiff');

close(f)
end
