% ROItype can be 'BigROI' or 'SmallROI'

function CombinedCorrMatrix(dataname, Acquisition, SaveDir, ROItype)

%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end
savename = ['timecourses_' dataname];

if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

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
    allrois =  {'VisualROI_R', 'SensoryROI_R', 'AuditoryROI_R', 'MotorROI_R',...
        'RetrosplenialROI_R', 'VisualROI_L','SensoryROI_L', 'AuditoryROI_L',...
        'MotorROI_L','RetrosplenialROI_L' };
else
    disp('ROItype not recognized.');
    return
end

%% Make groups for cacl/nacl
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat');
groups = {'CaCl', 'NaCl'};
RecordingOverview.Group = categorical(RecordingOverview.Group);
idx = RecordingOverview.Group == 'CaCl';

%% Start going per group, per mouse
for index = 1:size(groups,2) % Go per group, nacl/cacl
    group = groups{index};
    if matches(group, 'CaCl')
        Mousegroup = RecordingOverview(idx,:);
%         disp('cacl group')
    elseif matches(group, 'NaCl')
        Mousegroup = RecordingOverview(~idx, :);
%         disp('nacl group')
    else
        disp('something is going wrong')
        return
    end
    
    allcorrmatrices = zeros(size(Mousegroup,1), size(allrois,2), size(allrois,2));   %mice   x   roixroi (corrmatrix)
    
    for ind = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{ind};
        eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
        SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
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
        else
            Timecourses = NaN(10, 9000);
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
    
    saveas(gcf, [SaveDir '/CombinedCorrMatrices/' dataname '_' Acquisition '_' group '.tiff']);
    saveas(gcf, [SaveDir '/CombinedCorrMatrices/' dataname '_' Acquisition '_' group '.eps'], 'epsc');
    
    close(f)
    
    eval([group '_allcorrmatrices = allcorrmatrices;']);
    
end

%% plot Differences between groups
% z transform fisher
NaCl_allcorrmatrices = atanh(NaCl_allcorrmatrices);
CaCl_allcorrmatrices = atanh(CaCl_allcorrmatrices);

diff = mean(NaCl_allcorrmatrices, 1, 'omitnan') - mean(CaCl_allcorrmatrices, 1, 'omitnan');
diff = reshape(diff, size(allrois,2),size(allrois,2));

%% plot
f = figure('InvertHardcopy','off','Color',[1 1 1]);
ax = gca;
data = tril(diff);
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
load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
colormap(NL)
colorbar
f.Position = [10 10 1500 1500]; %for size of screen before saving

title(['Difference NaCl - CaCl ' Acquisition])

saveas(gcf, [SaveDir '/CombinedCorrMatrices/' dataname '_' Acquisition '_Difference.tiff']);
saveas(gcf, [SaveDir '/CombinedCorrMatrices/' dataname '_' Acquisition '_Difference.eps'], 'epsc');

close(f)

end




