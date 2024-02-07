% %% Check the spread of the correlation of a seed

% This is not the best way of doing it. It would be better if you have a
% big table or structure with the seedspreads for each  mouse so you can
% easily select the mice you need.

function CombinedSeedSpread(dataname, Acquisition, Grouping, Overwrite)

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

if ~exist([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep], 'dir')
    mkdir([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep])
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');

%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

[groups, ~] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end

disp(['Combined Seed spread ' dataname ' ' Acquisition ' ' Grouping]);

%% make table
[overviewtable] = MakeTable(dataname, Overwrite);

%% Plot groups together
plotcolours = [[0, 0.4470, 0.7410], [0.929, 0.694, 0.125], [0.635, 0.078, 0.184], [0.494, 0.184, 0.556], ...
    [0.3010 0.7450 0.9330],[1, 0.9, 0.1] ,[0.9, 0.1, 0.1] , [0.75, 0, 0.75]];
nrofROI = 8;

% overviewtable.idx = grp2idx(overviewtable.ROI); %this gives the groups based on alphabet, so sort the ROI labels as well:
codesdone = [];
labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};
nrofrounds = 7; %how many rounds of calculated spread you want to plot

for ind1 = 1:size(groups,1)
    group1 = groups{ind1}; %determine group 1
    codegroup1 = str2double(groups(ind1+size(groups,1)));

    for ind2 = 1:size(groups,1)
        group2 = groups{ind2}; %determine group 2
        codegroup2 = str2double(groups(ind2+size(groups,1)));
        
        paircode = codegroup1 + codegroup2;
        
        %don't compare the group with itself, don't do nonsense comparison
        %(sham male vs cacl female) and dont do if already done.
        if matches(group1, group2) || paircode == 33 || any(codesdone == paircode)
            continue
        end

        %% group 1
        %get table 8 x 15 (roi x steps)
        eval(['grouptable = overviewtable(overviewtable.' Grouping ' == group1,:);'])
        grouptable = grouptable(grouptable.Acquisition == Acquisition,:);
        Ngroup1 = num2str(size(unique(grouptable.Mouse),1));
        
        toplot = nan(length(labels), 15);
        stderror = nan(length(labels), 15);
        SEM = nan(length(labels), 15);
        
        for indroi = 1:length(labels)
            indexroi = grouptable.ROI == labels{indroi};
            temptable = grouptable.CorrSpread(indexroi,:);
            toplot(indroi,:) = mean(temptable, 1, 'omitnan');
            stderror(indroi,:) = std(temptable, 0, 1, 'omitnan');
            SEM(indroi,:) = std(temptable, 0, 1, 'omitnan')/sqrt(size(temptable,1));            
        end
        toplot = toplot';
        SEM = SEM';
        stderror = stderror';
        clear indroi indexroi grouptable temptable 

        %start plot
        f = figure('InvertHardcopy','off','Color',[1 1 1]);
        hold on

        % plot group 1 per roi
        indcol = 1;
        for ind = 1:nrofROI
            %1:7 is the number of spreadcircles you want to do
%             % no error bars:
%             plot(toplot(1:nrofrounds, ind),'LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
%             %std:
%             errorbar(toplot(1:nrofrounds, ind), stderror(1:nrofrounds, ind), 'LineWidth',2, 'Color', plotcolours(indcol:indcol+2))
            %SEM:
            errorbar(toplot(1:nrofrounds, ind), SEM(1:nrofrounds, ind), 'LineWidth',2, 'Color', plotcolours(indcol:indcol+2))
            indcol = indcol+3;
        end
        clear toplot indcol ind

        %% group 2
        eval(['grouptable = overviewtable(overviewtable.' Grouping ' == group2,:);'])
        grouptable = grouptable(grouptable.Acquisition == Acquisition,:);
        Ngroup2 = num2str(size(unique(grouptable.Mouse),1));
        toplot = nan(length(labels), 15);
        stderror = nan(length(labels), 15);
        SEM = nan(length(labels), 15);
        
        for indroi = 1:length(labels)
            indexroi = grouptable.ROI == labels{indroi};
            temptable = grouptable.CorrSpread(indexroi,:);
            toplot(indroi,:) = mean(temptable, 1, 'omitnan');
            stderror(indroi,:) = std(temptable, 0, 1, 'omitnan');
            SEM(indroi,:) = std(temptable, 0, 1, 'omitnan')/sqrt(size(temptable,1));            
        end
        toplot = toplot';
        SEM = SEM';
        stderror = stderror';
        clear indroi indexroi grouptable temptable 

        % plot group 2 per roi
        indcol = 1;
        for ind = 1:nrofROI
%             % no error bars
%             plot(toplot(1:nrofrounds, ind),'--','LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
%             %std:
%             errorbar(toplot(1:nrofrounds, ind), stderror(1:nrofrounds, ind), '--', 'LineWidth',2, 'Color', plotcolours(indcol:indcol+2))
            %SEM:
            errorbar(toplot(1:nrofrounds, ind), SEM(1:nrofrounds, ind), '--', 'LineWidth',2, 'Color', plotcolours(indcol:indcol+2))
            indcol = indcol+3;
        end
        clear toplot indcol ind

        %% Make pretty
        axes1 = gca;
        hold(axes1,'on');
        set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
        legend([labels labels], 'Location', 'southwest', 'NumColumns', 2)
        % f.Position = [10 10 1000 1000]; %for size of screen before saving
        f.Position = [10 10 1500 1500];
        ylim([0.5 1])
        xlabel('Radius')
        xticklabels({'Seed (3)', 'Round 1 (4-13)', 'Round 2 (14-23)','Round 3 (24-33)',...
            'Round 4 (34-43)','Round 5 (44-53)','Round 6 (54-63)'})
        xtickangle(45)
        title(['Spread of correlation ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none')
        subtitletext = {[group1 ' — vs ' group2 ' --'] ,...
            ['n = ' Ngroup1 ' vs n = ' Ngroup2], ...
            ['Error bars = SEM']};
        subtitle(subtitletext);
        
        % also saves when overwrite is 0, not ideal
        saveas(gcf, [SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group1 '-' group2 '-' Acquisition  '.tiff'], 'tiff');
        saveas(gcf, [SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group1 '-' group2 '-' Acquisition  '.eps'], 'epsc');
    
        close(f)
        codesdone = [codesdone, paircode];
        
    end
end

end % of function




function [overviewtable] = MakeTable(dataname, Overwrite)

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

% Go with both grouping variables, so you can combine them later if need be
Grouping = {'Group', 'Sex'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

Acquisitions = {'A1', 'A2', 'A3'};

%% Make table with overview
% make starting table
varTypes = {'cell', 'categorical', 'categorical', 'categorical', 'categorical', ...
    'categorical', 'single'};
varNames = {'Mouse','Acquisition','ROI','Combi','Group','Sex','CorrSpread'};
overviewtable = table('Size', [1 size(varNames,2)], 'VariableTypes', varTypes, 'VariableNames', varNames);
overviewtable.Mouse = 'dummy';
overviewtable.Acquisition = 'A1';

labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};

%check if table already exists. If so, and you dont have Overwrite on 1,
%load the table, go into the code anyway, jus skip if you already did the
%mouse and the acquisition
if exist([SaveDir '/SeedSpread/SeedSpreadTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 0
    load([SaveDir '/SeedSpread/SeedSpreadTable_' dataname(1:end-4) '.mat'], 'overviewtable');
elseif exist([SaveDir '/SeedSpread/SeedSpreadTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 1
    disp('Seed spread table already done, OVERWRITING MAT FILES')
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
            
            if matches(Mouse, 'M23') && matches(dataname, 'hemoCorr_fluo.dat')
                continue            
            elseif sum(temp) && sum(overviewtable(temp,:).Acquisition == Acquisition)
%                 disp(['Mouse ' Mouse ' ' Acquisition ' already done - skipped'])
                continue
            end
            
            eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
            SaveFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
            
            [CorrSpread] = SingleSubjectSeedSpread(SaveFolder, dataname);
            
            % take out auditory
            CorrSpread(3,:) = [];
            CorrSpread(7,:) = [];
            
            %build table single mouse
            tablemouse = table;
            tablemouse.Mouse = cellstr(repmat(Mousegroup.Mouse{indmouse}, size(labels,2),1));
            tablemouse.Acquisition = categorical(cellstr(repmat(Acquisition, size(labels,2),1)));
            tablemouse.ROI = categorical(labels');
            tablemouse.Combi = categorical(cellstr(repmat(group, size(labels, 2),1)));
            tablemouse.Group = categorical(cellstr(repmat(Mousegroup.Group(indmouse), size(labels, 2),1)));
            tablemouse.Sex = categorical(cellstr(repmat(Mousegroup.Sex(indmouse), size(labels, 2),1)));
            tablemouse.CorrSpread = CorrSpread;
            
            %add table to general table
            overviewtable = [overviewtable; tablemouse];
            
            clear tablemouse ConnValuesMouse PercentagesMouse DataFolder SaveFolder idx nrpixels
        end % of mice    
    end % of group
end % of acq

temp = matches(overviewtable.Mouse, 'dummy');
if sum(temp)
    overviewtable(temp,:) = [];
end

save([SaveDir '/SeedSpread/SeedSpreadTable_' dataname(1:end-4) '.mat'], 'overviewtable');

end

