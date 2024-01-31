%Method is 'std' or 'mean'
%Grouping is {'Sex'}, {'Group'}, {'Group','Sex'} or {'Combi'}
function CombinedGCaMPFluctuations(Acquisition, dataname, Grouping, Overwrite)

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

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
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

disp(['Combined GCaMP fluctuations ' dataname ' ' Acquisition ' ' Grouping]);

%% get table
[overviewtable] = MakeTable(dataname, Overwrite);
    
%% Plot boxplot
%% all 3 acq
[f, t] = MakeBoxplot123(Grouping, groups, overviewtable, [0 0.025], 'ROI', 'StdAct');

xlabel(t, 'Region of Interest', 'interpreter', 'none','FontSize',15,'FontWeight','bold')
ylabel(t, 'Fluctuation (std)', 'interpreter', 'none','FontSize',15,'FontWeight','bold');
title(t, ['Fluctuation ' dataname(1:end-4)], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
f.Position = [10 10 1500 1000]; %for size of screen before saving


if ~exist([SaveDir '/Fluctuations/Combined/' dataname(1:end-4)], 'dir')
    mkdir([SaveDir '/Fluctuations/Combined/' dataname(1:end-4)])
end

saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname(1:end-4) filesep 'boxplot_' Grouping '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname(1:end-4) filesep 'boxplot_' Grouping '.eps'], 'epsc');

close(f)

%% per acquisition
% get right Acquisition
tempacqindex = overviewtable.Acquisition == Acquisition;
overviewtable = overviewtable(tempacqindex,:);

[f, t] = MakeBoxplot(Grouping, groups, overviewtable, [0 0.025], 'ROI', 'StdAct');

if matches(Grouping, 'Sex')
    xlabel(t, 'Region of Interest', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(t, 'Fluctuation (std)', 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(t, ['Fluctuation ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
else
    xlabel('Region of Interest', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel('Fluctuation (std)', 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(['Fluctuation ' dataname(1:end-4) ' ' Acquisition], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
end

f.Position = [10 10 1500 1000]; %for size of screen before saving


% if ~exist([SaveDir '/Fluctuations/Combined/' dataname(1:end-4)], 'dir')
%     mkdir([SaveDir '/Fluctuations/Combined/' dataname(1:end-4)])
% end

saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname(1:end-4) filesep 'boxplot_' Acquisition '_' Grouping '.eps'], 'epsc');

close(f)

end




function [overviewtable] = MakeTable(dataname, Overwrite)

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
SaveDir = '/media/mbakker/GDrive/P2/GCaMP';

% Go with both grouping variables, so you can combine them later if need be
Grouping = {'Group', 'Sex'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

% Grouping = 'Combi';
Acquisitions = {'A1', 'A2', 'A3'};

%% Make table with overview
% make starting table
varTypes = {'cell', 'categorical', 'categorical', 'categorical', 'categorical', ...
    'categorical', 'single'};
varNames = {'Mouse','Acquisition','ROI','Combi','Group','Sex','StdAct'};
overviewtable = table('Size', [1 size(varNames,2)], 'VariableTypes', varTypes, 'VariableNames', varNames);
overviewtable.Mouse = 'dummy';
overviewtable.Acquisition = 'A1';

labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};

%check if table already exists. If so, and you dont have Overwrite on 1,
%load the table
if exist([SaveDir '/Fluctuations/FluctTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 0
    load([SaveDir '/Fluctuations/FluctTable_' dataname(1:end-4) '.mat'], 'overviewtable');
elseif exist([SaveDir '/Fluctuations/FluctTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 1
    disp('Fluctuations table already done, OVERWRITING MAT FILES')
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
            
            eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
            SaveFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
            %check if mouse is 23 (no fluo act), if there is a timecourse
            %and if the mouse is already in the table
            if matches(Mouse, 'M23') && matches(dataname, 'hemoCorr_fluo.dat')
                continue
            elseif ~exist([SaveFolder 'timecourses_' dataname(1:end-4) '_centroids.mat'], 'file')
                disp([Mouse ' does not have timecourses_' dataname '.mat. Run GetTimecourses first.'])                
                continue                 
            elseif sum(temp) && sum(overviewtable(temp,:).Acquisition == Acquisition)
%                 disp(['Mouse ' Mouse ' ' Acquisition ' already done - skipped'])
                continue
            end
            
            %get timecourses
            load([SaveFolder 'timecourses_' dataname(1:end-4) '_centroids.mat'], 'AllRois');
            Timecourses = cell2mat(AllRois(:,2));
            
            %get std of timecourses
            StdAct = movstd(Timecourses, 10, 0, 2, 'omitnan'); %0 is for w, default. [0 10] is to make the window forward moving
            
            if size(StdAct, 2) >= 9000
                StdAct = StdAct(:, 1:9000);
            else
                StdAct(:,end+1:9000) = missing;
            end
            
            %get rid of movement
            load([SaveFolder 'MovMask.mat'], 'MovMask');
            MovMask = MovMask(1:9000);
            StdAct = StdAct .* MovMask;
            StdAct(StdAct == 0) = nan;
            
            %get rid of auditory, take average 
            StdAct(3,:) = [];
            StdAct(7,:) = [];
            StdAct = mean(StdAct, 2, 'omitnan');
            
            %build table single mouse
            tablemouse = table;
            tablemouse.Mouse = cellstr(repmat(Mousegroup.Mouse{indmouse}, size(labels,2),1));
            tablemouse.Acquisition = categorical(cellstr(repmat(Acquisition, size(labels,2),1)));
            tablemouse.ROI = categorical(labels');
            tablemouse.Combi = categorical(cellstr(repmat(group, size(labels, 2),1)));
            tablemouse.Group = categorical(cellstr(repmat(Mousegroup.Group(indmouse), size(labels, 2),1)));
            tablemouse.Sex = categorical(cellstr(repmat(Mousegroup.Sex(indmouse), size(labels, 2),1)));
            tablemouse.StdAct = StdAct;
        
            %add mouse table to general table
            overviewtable = [overviewtable; tablemouse];
            
            clear tablemouse StdAct DataFolder SaveFolder idx Timecourses
        end % of mice    
    end % of group
end % of acq

temp = matches(overviewtable.Mouse, 'dummy');
if sum(temp)
    overviewtable(temp,:) = [];
end

save([SaveDir '/Fluctuations/FluctTable_' dataname(1:end-4) '.mat'], 'overviewtable');

end

