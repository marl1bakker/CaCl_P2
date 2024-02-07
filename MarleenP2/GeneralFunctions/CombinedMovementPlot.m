% imagingtype = {'GCaMP'} or {'Speckle'} or 'Both' or {'Both'}
% grouping = {'Combi'} or {'Sex'} or {'Group'}
% example: CombinedMovementPlot({'Combi'}, {'Speckle'}, 0);

function [statsresults] = CombinedMovementPlot(Grouping, imagingtype, Overwrite)

if ~exist('imagingtype', 'var')
    imagingtype = {'GCaMP'};
end
if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

%% Make 
f = figure('InvertHardcopy','off','Color',[1 1 1]);

if matches(imagingtype, 'Both')
    imagingtype = {'GCaMP', 'Speckle'};
    t = tiledlayout(1,2);
else
    t = tiledlayout("flow"); %if on 2023a or newer, do "horizontal", dont have to put it in the if/else then
end

for ind = 1:size(imagingtype,2)

    % get recordingoverview
if matches(imagingtype(ind), 'GCaMP')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
%     SaveDir = 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\MarleenP2\bla';
%     load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\MarleenP2\RecordingOverview.mat', 'RecordingOverview');

elseif matches(imagingtype(ind), 'Speckle')
    SaveDir = '/media/mbakker/GDrive/P2/Speckle';
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat', 'RecordingOverview');
%     SaveDir = 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\MarleenP2\bla';
%     load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\MarleenP2\RecordingOverview_Speckle.mat', 'RecordingOverview');
else
    disp('Imagingtype not recognized')
    return
end

%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

if size(Grouping, 2)>1 || matches(Grouping, 'Combi')
    Grouping = {'Combi'};
end

[groups, ~] = GroupVariables(RecordingOverview, Grouping);
Grouping = char(Grouping);

disp(['Combined movement plotting ' Grouping]);

%% Load .mat files of fluctuations
% Start going per group, per mouse
[overviewtable] = MakeTable(imagingtype{ind}, RecordingOverview, SaveDir, Overwrite);
if matches(imagingtype(ind), 'GCaMP')
    nrofframes = 9000; %hardcoded
elseif matches(imagingtype(ind), 'Speckle')
    nrofframes = 6000; %hardcoded
end

%% Plot boxplot    
nexttile
[t] = MakeBoxplotMovement(Grouping, groups, overviewtable, [0 nrofframes/3], 'Acquisition', 'Movement');

temp = gca;
sub = {temp.Subtitle.String, 'Outliers:'};
if any(overviewtable.Movement > (nrofframes/3))
    disp('*****')
    disp('MORE THAN 1/3 MOVEMENT!!')
    outliers = overviewtable(overviewtable.Movement>(nrofframes/3),:);
    % disp(outliers)

    for index = 1:size(outliers,1)
        sub = [sub, {[outliers.Mouse{index} ' ' char(outliers.Acquisition(index)) ' ' ...
            char(outliers.Combi(index)) ', nr of frames moved: ' ...
            num2str(outliers.Movement(index))]}];
    end
else
    sub = [sub, {'No outliers'}];
end
subtitle(sub)

if matches(Grouping, 'Sex')
    xlabel(t, 'Acquisition', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(t, ['Nr of frames moved (out of ' num2str(nrofframes) ')'], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(t, ['Movement ' imagingtype{ind}], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
else
    xlabel('Acquisition', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(['Nr of frames moved (out of ' num2str(nrofframes) ')'], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(['Movement ' imagingtype{ind}], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
end

%% statistics
eval(['statsresults.' imagingtype{ind} '= MovementStats(overviewtable);']);

end

f.Position = [10 50 1200 800]; 
leg = legend({'CaCl Female', 'CaCl Male', 'Sham Female', 'Sham Male', ''}, 'NumColumns', 4);
leg.Layout.Tile = 'south';

if sum(matches(imagingtype, {'GCaMP', 'Speckle'})) == 2
    imagingtype = {'Both'};
end

saveas(gcf, [SaveDir filesep 'Movement' filesep 'boxplot_' Grouping '_' imagingtype{1} '.tiff'], 'tiff');
saveas(gcf, [SaveDir filesep 'Movement' filesep 'boxplot_' Grouping '_' imagingtype{1} '.eps'], 'epsc');

close(f)

end



%%
%%




function [overviewtable] = MakeTable(imagingtype, RecordingOverview, SaveDir, Overwrite)

if matches(imagingtype, 'GCaMP')
    nrofframes = 9000; %hardcoded
elseif matches(imagingtype, 'Speckle')
    nrofframes = 6000; %hardcoded
end

% Go with both grouping variables, so you can combine them later if need be
Grouping = {'Group', 'Sex'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);
Acquisitions = {'A1', 'A2', 'A3'};

%% Make table with overview
% make starting table
varTypes = {'cell', 'categorical', 'categorical', 'categorical', ...
    'categorical', 'double'};
varNames = {'Mouse','Acquisition','Combi','Group','Sex','Movement'};
overviewtable = table('Size', [1 size(varNames,2)], 'VariableTypes', varTypes, 'VariableNames', varNames);
overviewtable.Mouse = 'dummy';
overviewtable.Acquisition = 'A1';

labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};

if exist([SaveDir filesep 'Movement' filesep 'MovementTable_' imagingtype '.mat'], 'file') && Overwrite == 0
    load([SaveDir filesep 'Movement' filesep 'MovementTable_' imagingtype '.mat'], 'overviewtable');
elseif exist([SaveDir filesep 'Movement' filesep 'MovementTable_' imagingtype '.mat'], 'file') && Overwrite == 1
    disp('Movement table already done, OVERWRITING MAT FILES')
end

for indacq = 1:size(Acquisitions, 2)
    Acquisition = Acquisitions{indacq};
    
    for indgroup = 1:size(groups,1) % Go per group
        group = groups{indgroup};
        %         disp(group);
        eval(['idx = RecordingOverview.Combi == ''' group ''';'])
        Mousegroup = RecordingOverview(idx,:);
        
        for indmouse = 1:size(Mousegroup, 1) %go per mouse
            Mouse = Mousegroup.Mouse{indmouse};
            temp = matches(overviewtable.Mouse, Mouse);
            
            if sum(temp) && sum(overviewtable(temp,:).Acquisition == Acquisition)
                %                 disp(['Mouse ' Mouse ' ' Acquisition ' already done - skipped'])
                continue
            elseif matches(Mousegroup.SaveDirectory{indmouse}, 'empty')
                continue
            end
            
            eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
            indices = strfind(DataFolder, '-');
            SaveFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep DataFolder(indices(end-1)+1:end) 'CtxImg' filesep];
            if ~exist(SaveFolder, 'dir') %this is because speckle doesnt have an CtxImg folder, while gcamp does
                SaveFolder = SaveFolder(1:end-7);
            end
            
            if exist([SaveFolder 'MovMask.mat'], 'file')
                load([SaveFolder 'MovMask.mat'], 'MovMask');
            else
                disp([SaveFolder ' is missing MovMask. Run Preprocessing pipeline first.'])
                continue
            end
            
            if size(MovMask, 2) >= nrofframes
                MovMask = MovMask(:, 1:nrofframes);
            else
                MovMask(:,end+1:nrofframes) = missing;
            end
            
            %build table single mouse
            tablemouse = table;
            tablemouse.Mouse = cellstr(Mousegroup.Mouse{indmouse});
            tablemouse.Acquisition = categorical(cellstr(Acquisition));
            tablemouse.Combi = categorical(cellstr(group));
            tablemouse.Group = categorical(Mousegroup.Group(indmouse));
            tablemouse.Sex = categorical(Mousegroup.Sex(indmouse));
            tablemouse.Movement = sum(MovMask == 0);

            %add table to general table
            overviewtable = [overviewtable; tablemouse];
            
            clear tablemouse Movement MovMask DataFolder SaveFolder idx indices
        end % of mice
    end % of group
end % of acq

temp = matches(overviewtable.Mouse, 'dummy');
if sum(temp)
    overviewtable(temp,:) = [];
end

if ~exist([SaveDir filesep 'Movement'], 'dir')
    mkdir([SaveDir filesep 'Movement'])
end
save([SaveDir filesep 'Movement' filesep 'MovementTable_' imagingtype '.mat'], 'overviewtable');

end






%%
%%

function [t] = MakeBoxplotMovement(Grouping, groups, overviewtable, ylimvalues, x, y)

if ~exist('x','var')
    x = 'ROI';
end

%name the columns in your table that you want to use x and y
temp = find(matches(overviewtable.Properties.VariableNames, x));
overviewtable.Properties.VariableNames{temp} = 'x';

if exist('y', 'var')
    temp = find(matches(overviewtable.Properties.VariableNames, y));
    overviewtable.Properties.VariableNames{temp} = 'y';
end

overviewtable.idx = grp2idx(overviewtable.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(overviewtable.x))';

if matches(Grouping, 'Sex') % if you're grouping by sex

    % BOXPLOT
    t = tiledlayout(2,1);
    conditions = {'CaCl','Sham'};

    %go per tile, first CaCl then Sham
    for indcondition = 1:size(conditions,2)
        Cond = overviewtable(overviewtable.Group == conditions{indcondition},:);

        axes = nexttile(t);
        b = boxchart(Cond.idx, Cond.y, 'GroupByColor', Cond.Sex,...
            'LineWidth', 2, 'MarkerStyle', 'none');
        hold on

        if indcondition == 1
            b(1).SeriesIndex = 7;
            b(2).SeriesIndex = 1;
        elseif indcondition == 2
            b(1).SeriesIndex = 2;
            b(2).SeriesIndex = 6;
        end

        % SCATTER
        xaxisstep = 1/size(groups,1);
        xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;

        for indroi = 1:length(labels)
            currentROI = Cond.x == labels{indroi};

            for indgroup = 1:size(groups,1)
                currentgroup = Cond.Sex == groups{indgroup};
                currentindex = currentgroup.*currentROI;
                Cond.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
            end

            xaxisplacement = xaxisplacement + 1;
        end
        clear indroi indgroup currentgroup currentROI currentindex xaxisplacement

        hs = scatter(Cond.idx, Cond.y, 70, 'filled', 'jitter','on','JitterAmount',0.02);

        % MAKE PRETTY
        ylim(ylimvalues);
        set(axes,'FontSize',20,'FontWeight','bold','LineWidth',2);
        legend('Female','Male','', 'Location', 'northeast', 'NumColumns', 2)
        hs.MarkerFaceColor = [0 0 0];
        hs.MarkerFaceAlpha = 0.3;
        hs.Marker;
        xticks(1:length(labels));
        xticklabels(labels)
        xlim([0.2 length(labels)+0.7])
        title(conditions{indcondition}, 'FontSize', 18)

        subtitleN = [];
        for indgroup = 1:size(groups, 1)
            ngroup = sum(Cond.Sex == groups{indgroup});
            ngroup = ngroup/length(labels);
            subtitleN = [subtitleN groups{indgroup} ' N = ' num2str(ngroup) ' -- '];
        end
        subtitle(axes, subtitleN(1:end-4), 'FontSize', 15)
    end

else

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
    % legend(b, 'Location', 'northeast', 'NumColumns', 2);
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
    subtitle(subtitleN(1:end-4), 'FontSize', 10)
    t = b; %to give return
end
end
