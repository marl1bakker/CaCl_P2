%Method is 'std' or 'mean'
function CombinedGCaMPFluctuations(Acquisition, dataname, Overwrite, Method, SaveDir)

if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

disp(['Combined Seed Spread ' Acquisition])

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end
%
% if ~exist('GSR', 'var')
%     GSR = 0;
% end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist('Method', 'var')
    Method = 'std';
end


% Make groups for cacl/nacl
groups = {'CaCl', 'NaCl'};
load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
RecordingOverview.Group = categorical(RecordingOverview.Group);
indcacl = RecordingOverview.Group == 'CaCl';

% Go per group, nacl/cacl
for indgroup = 1:size(groups,2)
    group = groups{indgroup};
    
    if exist([SaveDir '/Fluctuations/Combined/' group dataname(1:end-4) '_' Acquisition '.mat'], 'file') && ...
            Overwrite == 0
        disp([group ' already done, will load saved file.'])
        load([SaveDir '/Fluctuations/Combined/' group dataname(1:end-4) '_' Acquisition '.mat'], ['AllFlucts' group])
        savefile = 0; %dont save again if you already loaded it
        eval([group 'N = size(AllFlucts' group ', 3);'])
        continue
    elseif exist([SaveDir '/Fluctuations/Combined/' group dataname(1:end-4) '_' Acquisition '.mat'], 'file') && ...
            Overwrite == 1
        disp([group ' already done, WILL OVERWRITE.'])
    end
    savefile = 1;
    
    if matches(group, 'CaCl')
        Mousegroup = RecordingOverview(indcacl,:);
        disp('cacl group')
        CaClN = 0;
    elseif matches(group, 'NaCl')
        Mousegroup = RecordingOverview;
        Mousegroup(indcacl, :) = [];
        disp('nacl group')
        NaClN = 0;
    else
        disp('something is going wrong, group not recognized')
        return
    end

  
    AllFlucts = NaN(10, 9000, size(Mousegroup, 1)); % seeds, Timecourses(mask, timecourse, name), mice 
    
    for indmouse = 1:size(Mousegroup, 1) %go per mouse
        Mouse = Mousegroup.Mouse{indmouse};
        disp(Mouse)
        eval(['DataFolder = [Mousegroup.' Acquisition '{indmouse} filesep];']);
        SaveFolder = [SaveDir filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
        
        if exist([SaveFolder 'timecourses_' dataname '.mat'], 'file')
            load([SaveFolder 'timecourses_' dataname '.mat'], 'AllRois');
        else
            disp([Mouse ' does not have timecourses_' dataname '.mat. Run GetTimecourses first.'])
            AllFlucts(:, :, indmouse) = NaN(10,9000);

            continue
        end
        
        Timecourses = cell2mat(AllRois(:,2));
        
        if matches(Method, 'std')
            StdAct = movstd(Timecourses, 10, 0, 2, 'omitnan'); %0 is for w, default. [0 10] is to make the window forward moving
            %             elseif matches(Method, 'mean')
            % %                 StdAct = movmean(Timecourses, 10, 2, 'omitnan');
            %                 StdAct = Timecourses;
        end %if chosen mean method, the names are a bit weird since they keep referring to std, but the method is valid
        
        if size(StdAct, 2) >= 9000
            StdAct = StdAct(:, 1:9000);
        else
            StdAct(:,end+1:9000) = missing;
        end
        
        AllFlucts(:, :, indmouse) = StdAct;
        eval([group 'N = ' group 'N+1;'])
        
    end % of mice
    
    eval(['AllFlucts' group ' = AllFlucts;']);
    save([SaveDir '/Fluctuations/Combined/' group dataname(1:end-4) '_' Acquisition '.mat'], ...
        ['AllFlucts' group]); %save the matrix
    
end % of groups

%% plot
nrofROI = 8; %HARDCODED

nacl = mean(AllFluctsNaCl, 3, 'omitnan');
nacl = reshape(nacl, size(AllFluctsNaCl, 1), size(AllFluctsNaCl, 2));
nacl = nacl';
nacl(:,3) = []; %take out auditory
nacl(:,7) = [];

cacl = mean(AllFluctsCaCl, 3, 'omitnan');
cacl = reshape(cacl, size(AllFluctsCaCl, 1), size(AllFluctsCaCl, 2));
cacl = cacl';
cacl(:,3) = [];
cacl(:,7) = [];


plotcolours = [[0, 0.4470, 0.7410], [0.929, 0.694, 0.125], [0.635, 0.078, 0.184], [0.494, 0.184, 0.556], ...
    [0.3010 0.7450 0.9330],[1, 0.9, 0.1] ,[0.9, 0.1, 0.1] , [0.75, 0, 0.75]];
x = linspace(0,10,9000); %0 tot 10 min, 9000 frames
f = figure('InvertHardcopy','off','Color',[1 1 1]);
hold on

% plot NaCl
indcol = 1;
for ind = 1:nrofROI
    y = nacl(:,ind);
    y = movmedian(y, 1000, 1);
    plot(x,y,'LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); 
    indcol = indcol+3;
end

% plot CaCl
indcol = 1;
for ind = 1:nrofROI
    y = cacl(:,ind);
    y = movmedian(y, 1000, 1);
    plot(x,y,'--','LineWidth',2, 'Color', plotcolours(indcol:indcol+2)); 
    indcol = indcol+3;
end

% make pretty
axes1 = gca;
hold(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);

f.Position = [10 10 1500 500]; %for size of screen before saving

labels = {'Vis R', 'Sen R', 'Mot R', 'Ret R', 'Vis L', 'Sen L', 'Mot L', 'Ret L'};
legend(labels)

title(['Standard deviation of activity ' dataname ' ' Acquisition], 'interpreter', 'none')
subtitle(['N NaCl = ' num2str(NaClN) ' --- N CaCl = ' num2str(CaClN)]);

%save
if savefile == 1
    saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname '_' Acquisition  '.tiff'], 'tiff');
    saveas(gcf, [SaveDir '/Fluctuations/Combined/' dataname '_' Acquisition  '.eps'], 'epsc');
end

close(f)

end

























