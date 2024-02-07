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


%% Get correlation matrixes of the groups
% Single group Corr matrix
% Start going per group, per mouse
for index = 1:size(groups,1) % Go per group
    
    group = groups{index};
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])  
    Mousegroup = RecordingOverview(idx,:);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    
    allcorrmatrices = NaN(size(Mousegroup,1), size(allrois,2), size(allrois,2));   %mice   x   roixroi (corrmatrix)
    
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
%             imagesc(CorrMat, [0 1]);
%             disp(Mouse)
        end
        
        allcorrmatrices(ind,:,:) = CorrMat;
        clear Timecourses CorrMat AllRois MovMask
        
    end
    
    allcorrmatrices = reshape(allcorrmatrices,[],64);
    allcorrmatrices(isnan(allcorrmatrices(:,1)),:) = []; % get rid of nan mouse
    AvCorrMat = mean(allcorrmatrices, 1, 'omitnan');
    allcorrmatrices = reshape(allcorrmatrices,[],8,8);

    % Plot per group
    f = PlotCorrMatrix(AvCorrMat, allrois, [0 1]);

    title([dataname ' ' group ' ' Acquisition], 'interpreter', 'none')
    eval(['ngroup = ' group 'N;'])
    subtitle(['n = ' num2str(ngroup)])
    
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname(1:end-4) filesep dataname(1:end-4) '_' Acquisition '_' group '.tiff']);
    saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/PerGroup/' dataname(1:end-4) filesep dataname(1:end-4) '_' Acquisition '_' group '.eps'], 'epsc');
    
    close(f)
    allcorrmatrices(allcorrmatrices == 1) = NaN;
    eval([group '_allcorrmatrices = allcorrmatrices;']);
end

clear allcorrmatrices AvCorrMat ax ay data f group indx ind index Mouse ngroup NL Mousegroup Possibilities



%% Get differences between groups, get stats on them
% z transform fisher
for ind = 1:size(groups, 1)
   eval([ groups{ind} '_allcorrmatrices = atanh(' groups{ind} '_allcorrmatrices);'])
   eval([ groups{ind} '_allcorrmatrices = reshape(' groups{ind} '_allcorrmatrices, [], 64);'])
end

diff = [];
groupscompared = {};
codesdone = [];
allpvalues = [];

for ind1 = 1:size(groups,1)
    group1 = groups{ind1}; %determine group 1
    codegroup1 = str2double(groups(ind1+size(groups,1)));
    pvalues = NaN(64,1);

    for ind2 = 1:size(groups,1)
        group2 = groups{ind2}; %determine group 2
        codegroup2 = str2double(groups(ind2+size(groups,1)));
        
        paircode = codegroup1 + codegroup2;
        
        %don't compare the group with itself, don't do nonsense comparison
        %(sham male vs cacl female) and dont do if already done.
        if matches(group1, group2) || paircode == 33 || any(codesdone == paircode)
            continue
        end
        
        for indROI = 1:64
            % non parametric,
            % adtest(CaCl_allcorrmatrices(:,1,3)), varying over ROI. This
            % particular one is not parametric, so ill do all of them
            % nonparametric. Cant do wilcoxon because that is paired, so will
            % do mann whitney u
            if eval(['sum(isnan(' group1 '_allcorrmatrices(:,indROI))) == size(' group1 '_allcorrmatrices, 1)']) || ...
                    eval(['sum(isnan(' group2 '_allcorrmatrices(:,indROI))) == size(' group2 '_allcorrmatrices, 1)'])
                continue
            end
            
            eval(['pvalues(indROI) = ranksum(' group1 '_allcorrmatrices(:,indROI), ' group2 '_allcorrmatrices(:,indROI));']);
            
        end
        
        eval(['diffadd = mean(' group1 '_allcorrmatrices, 1, ''omitnan'') - mean(' group2 '_allcorrmatrices, 1, ''omitnan'');']);
        diffadd = reshape(diffadd, size(allrois,2)*size(allrois,2),1);
        diff = [diff, diffadd];
        groupscompared = [groupscompared; groups(ind1), groups(ind2)];
        codesdone = [codesdone; paircode];
        allpvalues = [allpvalues, pvalues];
    end
end
clear group1 group2 codegroup1 codegroup2 ind idx ind1 ind2 indROI paircode diffadd codesdone pvalues 

%% Do FDR
allqvalues = NaN(size(allpvalues));
allpqcombi = NaN(8,8,size(allpvalues,1));

for indgroup = 1:size(allpvalues,2) %if you have more than 2 groups...
    %get only non-repeated pvalues
    pvalues = reshape(allpvalues(:,indgroup), 8,8);
    pvalues = tril(pvalues, -1);
    pvalues = reshape(pvalues, 64, 1);
    
    pvalues28 = [];
    for ind = 1:size(pvalues,1)
        if pvalues(ind) ~= 0
            pvalues28 = [pvalues28; pvalues(ind)];
        end
    end
    clear ind
    
    % FDR
    qvalues28 = mafdr(reshape(pvalues28, [], 1),'BHFDR', 'true');
    
    % get p and q for 64, without repeaters
    q = tril(ones(8), -1);
    p = tril(ones(8), -1);
    q = reshape(q, 64, 1);
    p = reshape(p, 64, 1);
    
    ind2 = 1;
    for ind = 1:64
        if q(ind) ~= 0
            q(ind) = qvalues28(ind2);
            p(ind) = pvalues28(ind2);
            %         disp('Samuel est le meilleur programeur au monde');
            ind2 = ind2 + 1;
        end
    end
    
    q = reshape(q,8,8);
    p = reshape(p,8,8);
    pqcombi = p' + q;
    allpqcombi(:,:,indgroup) = pqcombi;   
    
    q(q==0) = NaN;
    p(p==0) = NaN;       
    q = reshape(q,64,1);
    p = reshape(p,64,1);
    allqvalues(:,indgroup) = q;
    allpvalues(:,indgroup) = p;
end

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



%% Plot Connectivity graph of differences between groups

% for ind = 1:size(allqvalues,2)
%     
%     CorrelationNetworkGraph(qvalues, diff)
%     
% 
% end
pqcombi(pqcombi==0) = NaN;
if any(pqcombi<0.05)
    disp('SIGN VALUE, MAKE CODE FOR FUCKED UP GRAPH')
end

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

% function [f] = CorrelationNetworkGraph(statistics, difference)
% % set up
% load('/home/mbakker/P2_scripts/MarleenP2/LargerAreaAtlas.mat', 'Map')
% Map(isnan(Map)) = 11;
% 
% %% get centroids
% Centroids = [];
% 
% for ind = 1:10 %For every ROI, get centroid to plot
%     [X, Y] = meshgrid(1:size(Map,2), 1:size(Map,1));
%     
%     %Get mask from only the ROI
%     ROI = Map;
%     ROI(ROI == ind) = 100;
%     ROI(ROI < 100) = 0;
%     ROI(ROI == 100) = 1;
%     ROI(isnan(ROI)) = 0;
%     
%     iX = sum(X(:).*ROI(:))/sum(ROI(:));
%     iY = sum(Y(:).*ROI(:))/sum(ROI(:));
%     iX = round(iX);
%     iY = round(iY);
%     Centroids = [Centroids; iX, iY];
% % end
% 
% %To exclude Auditory cortex
% Centroids = [Centroids(1:2,:); Centroids(4:7,:); Centroids(9:10,:)];
% clear ind iX iY ROI X Y A
% 
% %% plot
% f = figure();
% imcontour(Map, 'Color', 'k')
% 
% end
