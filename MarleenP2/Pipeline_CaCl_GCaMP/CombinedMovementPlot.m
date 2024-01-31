% imagingtype = 'GCaMP' or 'Speckle'

function CombinedMovementPlot(Grouping, imagingtype, Overwrite)

if ~exist('imagingtype', 'var')
    imagingtype = 'GCaMP';
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if matches(imagingtype, 'GCaMP')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
elseif matches(imagingtype, 'Speckle')
    SaveDir = '/media/mbakker/GDrive/P2/Speckle';
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat', 'RecordingOverview');
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
[overviewtable] = MakeTable(imagingtype, Overwrite);
if matches(imagingtype, 'GCaMP')
    nrofframes = 9000; %hardcoded
elseif matches(imagingtype, 'Speckle')
    nrofframes = 6000; %hardcoded
end
%% Plot boxplot    
[f, t] = MakeBoxplot(Grouping, groups, overviewtable, [0 nrofframes/3], 'Acquisition', 'Movement');

if any(overviewtable.Movement > (nrofframes/3))
    disp('*****')
    disp('MORE THAN 1/3 MOVEMENT!!')
    outliers = overviewtable(overviewtable.Movement>(nrofframes/3),:);
    disp(outliers)
    temp = gca;
    sub = {temp.Subtitle.String, 'Outliers:'};
    
    for index = 1:size(outliers,1)
        sub = [sub, {[outliers.Mouse{index} ' ' char(outliers.Acquisition(index)) ' ' ...
            char(outliers.Combi(index)) ', nr of frames moved: ' ...
            num2str(outliers.Movement(index))]}];
    end
    subtitle(sub)
end

if matches(Grouping, 'Sex')
    xlabel(t, 'Acquisition', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(t, ['Nr of frames moved (out of ' num2str(nrofframes) ')'], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(t, ['Movement ' imagingtype], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
else
    xlabel('Acquisition', 'interpreter', 'none','FontSize',20,'FontWeight','bold')
    ylabel(['Nr of frames moved (out of ' num2str(nrofframes) ')'], 'interpreter', 'none','FontSize',20,'FontWeight','bold');
    title(['Movement ' imagingtype], 'interpreter', 'none','FontSize',20,'FontWeight','bold')
end

f.Position = [10 10 800 1000]; %for size of screen before saving

if ~exist([SaveDir '/Movement'], 'dir')
    mkdir([SaveDir '/Movement'])
end

saveas(gcf, [SaveDir '/Movement/boxplot_' Grouping '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/Movement/boxplot_' Grouping '.eps'], 'epsc');

close(f)

end




function [overviewtable] = MakeTable(imagingtype, Overwrite)

SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
if matches(imagingtype, 'GCaMP')
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
    nrofframes = 9000; %hardcoded
elseif matches(imagingtype, 'Speckle')
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_Speckle.mat', 'RecordingOverview');
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

if exist([SaveDir '/Movement/MovementTable_' imagingtype '.mat'], 'file') && Overwrite == 0
    load([SaveDir '/Movement/MovementTable_' imagingtype '.mat'], 'overviewtable');
elseif exist([SaveDir '/Movement/MovementTable_' imagingtype '.mat'], 'file') && Overwrite == 1
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
%             SaveFolder = [Mousegroup.SaveDirectory{indmouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
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

save([SaveDir '/Movement/MovementTable_' imagingtype '.mat'], 'overviewtable');

end

