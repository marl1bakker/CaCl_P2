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
savename = ['timecourses_' dataname];

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
            Timecourses = Timecourses(:,1:9000);
            Timecourses(3,:) = [];
            Timecourses(7,:) = [];
            
            eval([group 'N = ' group 'N+1;'])
        else
            Timecourses = NaN(8, 9000);
        end
        
        CorrMat = corr(Timecourses');
        allcorrmatrices(ind,:,:) = CorrMat;
        clear Timecourses CorrMat AllRois
        
    end
    
    AvCorrMat = mean(allcorrmatrices, 1, 'omitnan');
    
    %% Plot per group
    f = figure('InvertHardcopy','off','Color',[1 1 1]);
    ax = gca;
    data = tril(reshape(AvCorrMat, size(allrois,2), size(allrois,2)));
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN; 
    imagesc(ax, data, 'AlphaData', ~isnan(data), [0 1])
    ax.Box = 'off';
    axis image;
        
    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ay,'FontSize',15, 'FontWeight', 'bold', 'Linewidth', 2);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)
    load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
    colormap(NL)
    colorbar
    f.Position = [10 10 1500 1500]; %for size of screen before saving
    title([dataname ' ' group ' ' Acquisition], 'interpreter', 'none')
    eval(['ngroup = ' group 'N;'])
    subtitle(['n = ' num2str(ngroup)])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname '_' Acquisition '_' group '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname '_' Acquisition '_' group '.eps'], 'epsc');
    
    close(f)
    
    eval([group '_allcorrmatrices = allcorrmatrices;']);
%     eval(['Mice' group '= Mousegroup;']);
end

clear allcorrmatrices AvCorrMat ax ay data f group idx ind index Mouse ngroup NL

%% Differences between groups

% z transform fisher
for ind = 1:size(groups, 1)
   eval([ groups{ind} '_allcorrmatrices = atanh(' groups{ind} '_allcorrmatrices);'])
end

diff = [];
% titles = {};
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
%         titles = [titles; strcat(group1, '-', group2)];
        groupscompared = [groupscompared; groups(ind1), groups(ind2)];
        codesdone = [codesdone; paircode];
    end
end
clear group1 group2 codegroup1 codegroup2 ind1 ind2 paircode 

%% plot
for ind = 1:size(diff, 2) %how many comparisons do you have
%     groupscompared = titles{ind};
    data = diff(:,ind);
    data = reshape(data, size(allrois,2),size(allrois,2));
    
    f = figure('InvertHardcopy','off','Color',[1 1 1]);
    ax = gca;
    data = tril(data);
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN;
    imagesc(ax, data, 'AlphaData', ~isnan(data), [-0.5 0.5])
    ax.Box = 'off';
    axis image;
    
    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ay,'FontSize', 15, 'FontWeight', 'bold', 'Linewidth', 2);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)
    load('/home/mbakker/P2_scripts/MarleenP2/NL.mat');
    colormap(NL)
    colorbar
    f.Position = [10 10 1500 1500]; %for size of screen before saving
    
    t = ['Difference ' groupscompared{ind, 1} ' - ' groupscompared{ind,2}, ...
        ', ' Acquisition];
    title(t)
    eval(['N1 = num2str(' groupscompared{ind,1} 'N);']);
    eval(['N2 = num2str(' groupscompared{ind,2} 'N);']);
    subtitle(['n ' groupscompared{ind,1} ' = ' N1 ' --- n ' groupscompared{ind, 2} '= ' N2])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/Differences/' dataname '_' groupscompared{ind,1} '-' groupscompared{ind,2} '_' Acquisition '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/Differences/' dataname '_' groupscompared{ind,1} '-' groupscompared{ind,2} '_' Acquisition '.eps'], 'epsc');
    
    close(f)

end

%% stats



end %of fucntion

% %% Make a group column in RecordingOverview or use an existing one. 
% switch size(Grouping,2)
%   
%     case 0
%         disp('No grouping variable selected, function exited')
%         return
%         
%     case 1
%         groups = eval(['cellstr(unique(RecordingOverview.' Grouping{1} '));']);
%         
%     otherwise %more than 1 grouping variable - make a Combi group
%         RecordingOverview.Combi = cellstr(repmat(' ', size(RecordingOverview,1),1));
%         for ind = 1:size(Grouping, 2)
%             eval(['currentgroup = cellstr(RecordingOverview.' Grouping{ind} ');' ]);
%             RecordingOverview.Combi = strcat(RecordingOverview.Combi, currentgroup);
%         end
%         groups = unique(RecordingOverview.Combi);
%         Grouping = 'Combi';
%         RecordingOverview.Combi = categorical(RecordingOverview.Combi);
% end
% clear currentgroup ind indx Possibilities ROItype
