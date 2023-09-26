% %% Check the spread of the correlation of a seed
function CombinedSeedSpread(dataname, Acquisition, SaveDir, Overwrite)

% set up
if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end
%
% if ~exist('GSR', 'var')
%     GSR = 0;
% end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end


% Make groups for cacl/nacl
groups = {'CaCl', 'NaCl'};
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
RecordingOverview.Group = categorical(RecordingOverview.Group);
indcacl = RecordingOverview.Group == 'CaCl';

% Go per group, nacl/cacl
for indgroup = 1:size(groups,2)
    group = groups{indgroup};
    
    if exist([SaveDir '/SeedSpread/Combined/CorrSpreads' group dataname(1:end-4) '_' Acquisition '.mat'], 'file') && ...
            Overwrite == 0
        disp([group ' already done, will load saved file.'])
        continue
    elseif exist([SaveDir '/SeedSpread/Combined/CorrSpreads' group dataname(1:end-4) '_' Acquisition '.mat'], 'file') && ...
            Overwrite == 1
        disp([group ' already done, WILL OVERWRITE.'])
    end
    
    if matches(group, 'CaCl')
        Mousegroup = RecordingOverview(indcacl,:);
        disp('cacl group')
    elseif matches(group, 'NaCl')
        Mousegroup = RecordingOverview;
        Mousegroup(indcacl, :) = [];
        disp('nacl group')
    else
        disp('something is going wrong, group not recognized')
        return
    end
    
    AllCorrSpreads = NaN(10, size(Mousegroup, 1), 15); % seeds, mice, circle steps (fixed in SingleSubjectSeedSpread)
    
    for indmouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indmouse};
        disp(Mouse)
        eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
        SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
        CorrSpread = SingleSubjectSeedSpread(SaveFolder, dataname);
        AllCorrSpreads(:, indmouse, :) = CorrSpread;
        
    end % of mice
    
    eval(['CorrSpreads' group ' = AllCorrSpreads;']);
    save([SaveDir '/SeedSpread/Combined/CorrSpreads' group dataname(1:end-4) '_' Acquisition '.mat'], ...
        ['CorrSpreads' group]); %save the matrix
    
end % of groups



%% plot
%get data
if ~exist('CorrSpreadsCaCl', 'var')
    load([SaveDir '/SeedSpread/Combined/CorrSpreadsCaCl' dataname(1:end-4) '_' Acquisition '.mat'], ...
        'CorrSpreadsCaCl');
end
if ~exist('CorrSpreadsNaCl', 'var')
    load([SaveDir '/SeedSpread/Combined/CorrSpreadsNaCl' dataname(1:end-4) '_' Acquisition '.mat'], ...
        'CorrSpreadsNaCl');
end

% hardcoded for 8 ROI
nrofROI = 8;

nacl = mean(CorrSpreadsNaCl, 2, 'omitnan');
nacl = reshape(nacl, size(CorrSpreadsNaCl, 1), size(CorrSpreadsNaCl, 3));
nacl = nacl';
nacl(:,3) = []; %take out auditory
nacl(:,7) = [];

cacl = mean(CorrSpreadsCaCl, 2, 'omitnan');
cacl = reshape(cacl, size(CorrSpreadsCaCl, 1), size(CorrSpreadsCaCl, 3));
cacl = cacl';
cacl(:,3) = [];
cacl(:,7) = [];

plotcolours = [[0, 0.4470, 0.7410], [0.929, 0.694, 0.125], [0.635, 0.078, 0.184], [0.494, 0.184, 0.556], ...
        [0.3010 0.7450 0.9330],[1, 0.9, 0.1] ,[0.9, 0.1, 0.1] , [0.75, 0, 0.75]]; 

f = figure('InvertHardcopy','off','Color',[1 1 1]);
hold on

% plot NaCl
indcol = 1;
for ind = 1:nrofROI
    plot(nacl(1:6, ind),'LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
    indcol = indcol+3;
end

% plot CaCl
indcol = 1;
for ind = 1:nrofROI
    plot(cacl(1:6, ind),'--','LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); %only take the first six rounds of periphery
    indcol = indcol+3;
end

%make pretty
axes1 = gca;
hold(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
labels = {'Vis R', 'Sen R', 'Mot R', 'Ret R', 'Vis L', 'Sen L', 'Mot L', 'Ret L'};
legend([labels labels], 'Location', 'southwest', 'NumColumns', 2)
f.Position = [10 10 1000 1000]; %for size of screen before saving
ylim([0.5 1])
title(['Spread of correlation ' dataname(1:end-4) Acquisition])

saveas(gcf, [SaveDir '/SeedSpread/Combined/CorrSpread' dataname(1:end-4) '_' Acquisition  '.tiff'], 'tiff');
saveas(gcf, [SaveDir '/SeedSpread/Combined/CorrSpread' dataname(1:end-4) '_' Acquisition  '.eps'], 'epsc');

end % of function

