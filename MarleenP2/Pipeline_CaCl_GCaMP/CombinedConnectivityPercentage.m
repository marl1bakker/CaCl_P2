function CombinedConnectivityPercentage(dataname, Acquisition, Grouping, threshold, Overwrite)

%% set up
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

[groups, ~] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end

%% Make table with overview
[overviewtable] = MakeTable(dataname, threshold, Overwrite);

%% start plotting
overviewtable.Perc = overviewtable.Perc * 100;

% get right Acquisition
tempacqindex = overviewtable.Acquisition == Acquisition;
overviewtable = overviewtable(tempacqindex,:);

%which variable do you want to set as y:
temp = find(matches(overviewtable.Properties.VariableNames, 'Perc'));
overviewtable.Properties.VariableNames{temp} = 'y';

[f, t] = MakeBoxplot(Grouping, groups, overviewtable, [0 100]);

if matches(Grouping, 'Sex')
    xlabel(t, 'Region of Interest', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(t, ['% pixels with correlation > ' num2str(threshold)], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(t, ['Connectivity  ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
else
    xlabel('Region of Interest', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(['% pixels with correlation > ' num2str(threshold)], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(['Connectivity ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
end

f.Position = [10 10 1500 1000]; %for size of screen before saving

if ~exist([SaveDir '/ConnectivityPercentage/' dataname(1:end-4)], 'dir')
    mkdir([SaveDir '/ConnectivityPercentage/' dataname(1:end-4)])
end

saveas(gcf, [SaveDir '/ConnectivityPercentage/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/ConnectivityPercentage/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.eps'], 'epsc');

close(f)

end


function [overviewtable] = MakeTable(dataname, threshold, Overwrite)

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');

% Go with both grouping variables, so you can combine them later if need be
Grouping = {'Group', 'Sex'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);
Acquisitions = {'A1', 'A2', 'A3'};

%% Make table with overview
% make starting table
varTypes = {'cell', 'categorical', 'categorical', 'categorical', 'categorical', ...
    'categorical', 'double', 'double'};
varNames = {'Mouse','Acquisition','ROI','Combi','Group','Sex','Perc', 'Pixels'};
overviewtable = table('Size', [1 size(varNames,2)], 'VariableTypes', varTypes, 'VariableNames', varNames);
overviewtable.Mouse = 'dummy';
overviewtable.Acquisition = 'A1';

labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};

if exist([SaveDir '/ConnectivityPercentage/ConnPercTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 0
    load([SaveDir '/ConnectivityPercentage/ConnPercTable_' dataname(1:end-4) '.mat'], 'overviewtable');
elseif exist([SaveDir '/ConnectivityPercentage/ConnPercTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 1
    disp('Connectivity percentage table already done, OVERWRITING MAT FILES')
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
            
            % get connectivity and percentage values
            [ConnValuesMouse, nrpixels] = SingleSubjectConnectivityValues(SaveFolder, dataname, threshold);
            PercentagesMouse = ConnValuesMouse/nrpixels;
            
            ConnValuesMouse = ConnValuesMouse';
            ConnValuesMouse(3,:) = []; %take out auditory
            ConnValuesMouse(7,:) = []; %toplot is roi x mice
            
            PercentagesMouse = PercentagesMouse';
            PercentagesMouse(3,:) = [];
            PercentagesMouse(7,:) = [];
            
            %build table single mouse
            tablemouse = table;
            tablemouse.Mouse = cellstr(repmat(Mousegroup.Mouse{indmouse}, size(labels,2),1));
            tablemouse.Acquisition = categorical(cellstr(repmat(Acquisition, size(labels,2),1)));
            tablemouse.ROI = categorical(labels');
            tablemouse.Combi = categorical(cellstr(repmat(group, size(labels, 2),1)));
            tablemouse.Group = categorical(cellstr(repmat(Mousegroup.Group(indmouse), size(labels, 2),1)));
            tablemouse.Sex = categorical(cellstr(repmat(Mousegroup.Sex(indmouse), size(labels, 2),1)));
            tablemouse.Perc = PercentagesMouse;
            tablemouse.Pixels = ConnValuesMouse;
            
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

save([SaveDir '/ConnectivityPercentage/ConnPercTable_' dataname(1:end-4) '.mat'], 'overviewtable');

end