function CombinedCorrMatrix(SaveDir, Acquisition)

if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

groups = {'CaCl', 'NaCl'};

%get all possible roi
load('/home/mbakker/P2_scripts/Umit-release_Astrocyte-v1.5/GUI/DataViz/mouse_ctx_borders.mat')
allrois = [atlas.areatag atlas.area_longNames atlas.area_ClusterNames];
clear atlas
[~, I] = sort(allrois(1:24,3)); %sort Timecourses based on groups in allrois. do 24 instead of 48 so you keep left and right
I = [I;I+24];
        
for indexx = 1:size(groups,2) % Go per group, nacl/cacl
    group = groups{indexx};
    
    % Make groups for cacl/nacl
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat');
    RecordingOverview.Group = categorical(RecordingOverview.Group);
    idx = RecordingOverview.Group == 'CaCl';
    
    if matches(group, 'CaCl')
        Mousegroup = RecordingOverview(idx,:);
        disp('cacl group')
    elseif matches(group, 'NaCl')
        Mousegroup = RecordingOverview(~idx, :);
        disp('nacl group')
    else
        disp('something is going wrong')
        return
    end
    clear idx RecordingOverview
    
    allcorrmatrices = zeros(size(Mousegroup,1), size(allrois,1), size(allrois,1));   %mice   x   roixroi (corrmatrix)
    
    for ind = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{ind};
        eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
        SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        load([SaveFolder 'AllRoisFluo.mat']);
        
        % make sure everyone has the same roi
        roismouse = AllRoisFluo(:,3);
        missingrois = ismember(allrois(:,1), roismouse); %rois that are not present have a 0
        
        for index = find(~missingrois') %find places where missingrois is 0
            %         disp(index)
            insertrow = {logical(zeros(512, 512)), single(NaN(1, size(AllRoisFluo{1,2},2))), allrois(index)};
            AllRoisFluo = [AllRoisFluo(1:index-1,:); insertrow; AllRoisFluo(index:end,:)];
%             AllRoisHbO = [AllRoisHbO(1:index-1,:); insertrow; AllRoisHbO(index:end,:)];
        end
        clear DataFolder SaveFolder index insertrow missingrois roismouse Mouse
        
        Timecourses = cell2mat(AllRoisFluo(:,2));
        Timecourses = Timecourses(:,1:9000); %not really needed, but nice to have all same size
        Timecourses = Timecourses(I,:); %made I before loop, groups larger regions
        
        % all ROIs
        CorrMat = corr(Timecourses');
        allcorrmatrices(ind,:,:) = CorrMat;
        
        % clustered ROIs
        
        clear Timecourses CorrMat AllRoisFluo
        
    end
    
    AvCorrMat = mean(allcorrmatrices, 1, 'omitnan');
    
    % plot per group
    f = figure('InvertHardcopy','off','Color',[1 1 1]);
    ax = gca;
    data = tril(reshape(AvCorrMat,48,48));
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN; %if you want to get rid of the 1's in the diagonal row
    imagesc(ax, data, 'AlphaData', ~isnan(data), [-1 1])
    ax.Box = 'off';
    axis image;
    
    tags = allrois(:,1);
    tags = tags(I,:);
    
    yticks(1:48);
    yticklabels(tags);
    ay = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ay,'FontSize', 10, 'FontWeight', 'bold', 'Linewidth', 2);
    xticks(1:48);
    xticklabels(tags);
    xtickangle(90)
    load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
    colormap(NL)
    colorbar
    f.Position = [10 10 1500 1500]; %for size of screen before saving
    title([group ' ' Acquisition])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/' Acquisition '_' group '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/' Acquisition '_' group '.eps'], 'epsc');
    
    close(f)
    
    eval([group '_allcorrmatrices = allcorrmatrices;']);

end

%% plot Differences between groups

% z transform fisher
NaCl_allcorrmatrices = atanh(NaCl_allcorrmatrices);
CaCl_allcorrmatrices = atanh(CaCl_allcorrmatrices);

diff = mean(NaCl_allcorrmatrices, 1, 'omitnan') - mean(CaCl_allcorrmatrices, 1, 'omitnan');
diff = reshape(diff, 48,48);

f = figure('InvertHardcopy','off','Color',[1 1 1]);
ax = gca;
data = tril(diff);
data(data == 0 ) = NaN;
data(data == 1 ) = NaN; %if you want to get rid of the 1's in the diagonal row
imagesc(ax, data, 'AlphaData', ~isnan(data), [-1 1])
ax.Box = 'off';
axis image;

yticks(1:48);
yticklabels(tags);
ay = get(gca,'YTickLabel');
set(gca,'YTickLabel',ay,'FontSize', 10, 'FontWeight', 'bold', 'Linewidth', 2);
xticks(1:48);
xticklabels(tags);
xtickangle(90)
load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
colormap(NL)
colorbar
f.Position = [10 10 1500 1500]; %for size of screen before saving

title(['Difference NaCl - CaCl ' Acquisition])

saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/' Acquisition '_Difference.tiff']);
saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/' Acquisition '_Difference.eps'], 'epsc');
   
close(f)

%% pool regions together
clear AvCorrMat ax ay data diff f group groups ind indexx NL Mousegroup

allrois = allrois(I,:);
clusters = unique(allrois(:,3));

for ind = 1:size(clusters, 1)
    cluster = clusters(ind);
    index = find(strcmp(allrois(:,3), cluster));
    indexright = index(index<25);
    indexleft = index(index>24);
    
    mean(
    
end



end




