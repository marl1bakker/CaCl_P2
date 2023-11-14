% %% Check the spread of the correlation of a seed

% This is not the best way of doing it. It would be better if you have a
% big table or structure with the seedspreads for each  mouse so you can
% easily select the mice you need.

function CombinedSeedSpread(dataname, Acquisition, Grouping)

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

disp(['Combined Seed Spread ' Acquisition])

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

[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end
%% Single group Corr matrix - make or load data
% Start going per group, per mouse
for indgroup = 1:size(groups,1) % Go per group
    group = groups{indgroup};

    if exist([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group  '_' Acquisition '.mat'], 'file')
        load([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group  '_' Acquisition '.mat'], ['spreads' group]);
        % count the nans, if its 150 it is all nans so dont count the mouse
        eval(['nanspermouse = sum(sum(isnan(spreads' group '), 3),1);'])
        eval([group 'N = size(spreads' group ', 2) - length(find(nanspermouse == 150));']);
        clear nanspermouse
        continue
    end
    
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])  
    Mousegroup = RecordingOverview(idx,:);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    AllCorrSpreads = NaN(10, size(Mousegroup, 1), 15); % seeds, mice, circle steps (fixed in SingleSubjectSeedSpread)

    for indmouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indmouse};
        eval(['RawDataFolder = [Mousegroup.' Acquisition '{' num2str(indmouse) '} filesep];']);
        DataFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep RawDataFolder(end-5:end) 'CtxImg' filesep];
        
        CorrSpread = SingleSubjectSeedSpread(DataFolder, dataname);
        AllCorrSpreads(:, indmouse, :) = CorrSpread;    
        
        if ~isequal(sum(isnan(CorrSpread), 'all'), 10*15) %hardcoded 10 roi, 15 steps
            eval([group 'N = ' group 'N + 1;']);
        end
    end % of mice
    
    eval(['spreads' group ' = AllCorrSpreads;']);
    save([SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group  '_' Acquisition '.mat'], ...
        ['spreads' group]); %save the matrix
    
end % of groups
clear  AllCorrSpreads RawDataFolder DataFolder MouseGroup indmouse group

%% Plot per group
% for indgroup = 1:size(groups,1) 
%     if exist([SaveDir '/SeedSpread/Combined/' group '_' dataname(1:end-4) '_' Acquisition  '.eps'], 'epsc');
%         continue
%     end
%     
% end

%% Plot groups together
plotcolours = [[0, 0.4470, 0.7410], [0.929, 0.694, 0.125], [0.635, 0.078, 0.184], [0.494, 0.184, 0.556], ...
    [0.3010 0.7450 0.9330],[1, 0.9, 0.1] ,[0.9, 0.1, 0.1] , [0.75, 0, 0.75]];
nrofROI = 8;

codesdone = [];

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

        %start plotting group 1
        f = figure('InvertHardcopy','off','Color',[1 1 1]);
        hold on
        eval(['toplot = mean(spreads' group1 ', 2, ''omitnan'');']);
        toplot = reshape(toplot, 10,15); %hardcoded on nr of roi and nr of steps
        toplot = toplot';
        toplot(:,3) = []; %take out auditory
        toplot(:,7) = [];
        
        % plot group 1 per roi
        indcol = 1;
        for ind = 1:nrofROI
            plot(toplot(1:7, ind),'LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
            indcol = indcol+3;
        end
        clear toplot indcol ind
        
        %start plotting group 2
        eval(['toplot = mean(spreads' group2 ', 2, ''omitnan'');']);
        toplot = reshape(toplot, 10,15);
        toplot = toplot';
        toplot(:,3) = [];
        toplot(:,7) = [];
        
        % plot group 2 per roi
        indcol = 1;
        for ind = 1:nrofROI
            plot(toplot(1:7, ind),'--','LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
            indcol = indcol+3;
        end
        
        %Make pretty
        axes1 = gca;
        hold(axes1,'on');
        set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
        labels = {'Vis R', 'Sen R', 'Mot R', 'Ret R', 'Vis L', 'Sen L', 'Mot L', 'Ret L'};
        legend([labels labels], 'Location', 'southwest', 'NumColumns', 2)
        f.Position = [10 10 1000 1000]; %for size of screen before saving
        ylim([0.5 1])
        xlabel('Radius')
        xticklabels({'Seed (3)', 'Round 1 (4-13)', 'Round 2 (14-23)','Round 3 (24-33)',...
            'Round 4 (34-43)','Round 5 (44-53)','Round 6 (54-63)'})
        xtickangle(45)
        title(['Spread of correlation ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none')
        Ngroup1 = eval(['num2str(' group1 'N);']);
        Ngroup2 = eval(['num2str(' group2 'N);']);
        subtitletext = {[group1 ' — vs ' group2 ' --'] ,...
            ['n = ' Ngroup1 ' vs n = ' Ngroup2]};
        subtitle(subtitletext);
        
        % also saves when overwrite is 0, not ideal
        saveas(gcf, [SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group1 '-' group2 '-' Acquisition  '.tiff'], 'tiff');
        saveas(gcf, [SaveDir '/SeedSpread/Combined/' dataname(1:end-4) filesep group1 '-' group2 '-' Acquisition  '.eps'], 'epsc');
    
        close(f)
        codesdone = [codesdone, paircode];
        
    end
end

end % of function

