% ROItype can be 'BigROI' or 'SmallROI'
% Grouping can be 'Group', 'Sex', or 'Combi'

% This code is a bit crappy. If you select multiple groups for grouping, it
% also saves nonsense comparisons like caclmale with shamfemale. Haven't
% found a way yet to make this better. 

function CombinedCorrMatrix(dataname, Acquisition, Grouping, ROItype)

%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end
savename = ['timecourses_' dataname(1:end-4) '_centroids'];

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

disp(['Combined Correlation Matrix ' Acquisition])

if ~exist('ROItype', 'var')
    ROItype = 'BigROI';
end

% Choose right roitype
if matches(ROItype, 'SmallROI')
    %get all possible roi (small roi)
    % load('/home/mbakker/P2_scripts/Umit-release_Astrocyte-v1.5/GUI/DataViz/mouse_ctx_borders.mat')
    % allrois = [atlas.areatag atlas.area_longNames atlas.area_ClusterNames];
    % clear atlas
    % [~, I] = sort(allrois(1:24,3)); %sort Timecourses based on groups in allrois. do 24 instead of 48 so you keep left and right
    % I = [I;I+24];
elseif matches(ROItype, 'BigROI')
    allrois =  {'VisualROI_R', 'SensoryROI_R', 'MotorROI_R',...
        'RetrosplenialROI_R', 'VisualROI_L','SensoryROI_L', ...
        'MotorROI_L','RetrosplenialROI_L' };
else
    disp('ROItype not recognized.');
    return
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

%% Table


%% Single group Corr matrix
% Start going per group, per mouse
for index = 1:size(groups,1) % Go per group
    
    group = groups{index};
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])  
    Mousegroup = RecordingOverview(idx,:);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    
    allcorrmatrices = zeros(size(Mousegroup,1), size(allrois,2), size(allrois,2));   %mice   x   roixroi (corrmatrix)
    
    for ind = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{ind};
        eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
        SaveFolder = [Mousegroup.SaveDirectory{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
        if matches(Mouse, 'M23') && matches(dataname, 'hemoCorr_fluo.dat')
            continue %M23 has not enough fluorescence, so skip for fluo.
        end
        
        if exist([SaveFolder savename '.mat'], 'file')
            load([SaveFolder savename '.mat'], 'AllRois');
            
            % make sure everyone has the same roi (if going for small roi)
            %             roismouse = AllRoisFluo(:,3);
            %             missingrois = ismember(allrois(:,1), roismouse); %rois that are not present have a 0
            %
            %         for index = find(~missingrois') %find places where missingrois is 0
            %             %         disp(index)
            %             insertrow = {logical(zeros(512, 512)), single(NaN(1, size(AllRoisFluo{1,2},2))), allrois(index)};
            %             AllRoisFluo = [AllRoisFluo(1:index-1,:); insertrow; AllRoisFluo(index:end,:)];
            % %             AllRoisHbO = [AllRoisHbO(1:index-1,:); insertrow; AllRoisHbO(index:end,:)];
            %         end
            %         clear DataFolder SaveFolder index insertrow missingrois roismouse Mouse
            
            Timecourses = cell2mat(AllRois(:,2));
            load([SaveFolder 'MovMask.mat'], 'MovMask');
            Timecourses = Timecourses .* MovMask; %to remove movement. I checked, the difference between corr with and without movement is very small, mostly 0.01 differences. checked for one mouse (M15-A1)       
            Timecourses(Timecourses == 0) = nan;
            Timecourses = Timecourses(:,1:9000);
            Timecourses(3,:) = [];
            Timecourses(7,:) = [];

            eval([group 'N = ' group 'N+1;'])
        else
            Timecourses = NaN(8, 9000);
        end
        
        CorrMat = corr(Timecourses', 'rows','pairwise');
        
        if any(CorrMat < 0.4, 'all')
            imagesc(CorrMat, [0 1]);
            disp(Mouse)
        end
        
        allcorrmatrices(ind,:,:) = CorrMat;
        clear Timecourses CorrMat AllRois
        
    end
    
    AvCorrMat = mean(allcorrmatrices, 1, 'omitnan');
    
    %% Plot per group
    f = PlotCorrMatrix(AvCorrMat, allrois, [0 1]);

    title([dataname ' ' group ' ' Acquisition], 'interpreter', 'none')
    eval(['ngroup = ' group 'N;'])
    subtitle(['n = ' num2str(ngroup)])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname(1:end-4) filesep dataname(1:end-4) '_' Acquisition '_' group '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname(1:end-4) filesep dataname(1:end-4) '_' Acquisition '_' group '.eps'], 'epsc');
    
    close(f)
    
    eval([group '_allcorrmatrices = allcorrmatrices;']);
end

clear allcorrmatrices AvCorrMat ax ay data f group idx ind index Mouse ngroup NL

%% Get differences between groups
% z transform fisher
for ind = 1:size(groups, 1)
   eval([ groups{ind} '_allcorrmatrices = atanh(' groups{ind} '_allcorrmatrices);'])
end

diff = [];
groupscompared = {};
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
        
        eval(['diffadd = mean(' group1 '_allcorrmatrices, 1, ''omitnan'') - mean(' group2 '_allcorrmatrices, 1, ''omitnan'');']);
        diffadd = reshape(diffadd, size(allrois,2)*size(allrois,2),1);
        diff = [diff, diffadd];
        groupscompared = [groupscompared; groups(ind1), groups(ind2)];
        codesdone = [codesdone; paircode];
    end
end
clear group1 group2 codegroup1 codegroup2 ind1 ind2 paircode 

%% plot CorrMatrix of differences between groups
for ind = 1:size(diff, 2) %how many comparisons do you have

    f = PlotCorrMatrix(diff(:,ind), allrois, [-0.5 0.5]);
    
    t = ['Difference ' groupscompared{ind, 1} ' - ' groupscompared{ind,2}, ...
        ', ' Acquisition];
    title(t)
    eval(['N1 = num2str(' groupscompared{ind,1} 'N);']);
    eval(['N2 = num2str(' groupscompared{ind,2} 'N);']);
    subtitle(['n ' groupscompared{ind,1} ' = ' N1 ' --- n ' groupscompared{ind, 2} '= ' N2])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/Differences/' dataname(1:end-4) filesep dataname(1:end-4) '_' groupscompared{ind,1} '-' groupscompared{ind,2} '_' Acquisition '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/Differences/' dataname(1:end-4) filesep dataname(1:end-4) '_' groupscompared{ind,1} '-' groupscompared{ind,2} '_' Acquisition '.eps'], 'epsc');
    
    close(f)
end

%% Plot Connectivity graph

%% stats



end %of fucntion


function [f] = PlotCorrMatrix(data, allrois, axlims)
    f = figure('InvertHardcopy','off','Color',[1 1 1]);
    ax = gca;
    data = tril(reshape(data, size(allrois,2), size(allrois,2)));
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN; 
    imagesc(ax, data, 'AlphaData', ~isnan(data), axlims)
    ax.Box = 'off';
    axis image;
        
    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ay,'FontSize',15, 'FontWeight', 'bold', 'Linewidth', 2);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)
    load('/home/mbakker/P2_scripts/MarleenP2/NL.mat');
    colormap(NL)
    colorbar
    f.Position = [10 10 1500 1500]; %for size of screen before saving
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
varNames = {'Mouse','Acquisition','ROI','Combi','Group','Sex','CorrMatrix'};
overviewtable = table('Size', [1 size(varNames,2)], 'VariableTypes', varTypes, 'VariableNames', varNames);
overviewtable.Mouse = 'dummy';
overviewtable.Acquisition = 'A1';

labels = {'Vis-R', 'Sen-R', 'Mot-R', 'Ret-R', 'Vis-L', 'Sen-L', 'Mot-L', 'Ret-L'};

%check if table already exists. If so, and you dont have Overwrite on 1,
%load the table
if exist([SaveDir '/CombinedCorrMatrices/CorrMatTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 0
    load([SaveDir '/CombinedCorrMatrices/CorrMatTable_' dataname(1:end-4) '.mat'], 'overviewtable');
elseif exist([SaveDir '/CombinedCorrMatrices/CorrMatTable_' dataname(1:end-4) '.mat'], 'file') && Overwrite == 1
    disp('Correlation Matrix table already done, OVERWRITING MAT FILES')
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
