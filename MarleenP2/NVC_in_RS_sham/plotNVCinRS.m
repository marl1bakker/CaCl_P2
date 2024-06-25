% type can be:
% - 'normal', which takes hemoCorr_fluo and takes the average of
% each pixel, so that all pixels "weigh" the same.
% - 'nofilt', which takes the hemocorr_fluo data to find the activations,
% but then takes the non-filtered data to compute curves and specs
% - 'unweighted', same as normal but with a pixel with 100 activations
% weighing more than a pixel with 1. Takes all activations on a ROI,
% regardless of which pixel.

% ROIsavename can be:
% - 'LR', left vs right
% - 'BigROI', which takes all areas (visual, sensory etc.) except auditory
% - 'WholeBrain', no ROI, just the brain in general
%

function plotNVCinRS(Acq, type, ROIsavename)

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');
RecordingOverview = RecordingOverview(RecordingOverview.Group == 'Sham',:);

% note: outliersremoved and normal differs extremely little, but for the
% outliersremoved the dips and peaks are saved in the specs. 

if ~exist('type','var') || matches(type, 'normal')
    type = 'normal';
    datanamefluocurves = 'hemoCorr_fluo';
    NVCname = 'NVC_ROI';
elseif matches(type, 'OutliersRemoved')
    type = 'OutliersRemoved';
    datanamefluocurves = 'hemoCorr_fluo';
    NVCname = 'NVC_ROI_OL_removed';
elseif matches(type, 'nofilt')
    datanamefluocurves = 'fluo_nofilt';
    NVCname = ['NVC_ROI_' datanamefluocurves];
elseif matches(type, 'unweighted')
    NVCname = 'NVC_ROI_unweighted';
end

if ~exist('ROIsavename', 'var')
    ROIsavename = 'BigROI';
end

if matches(ROIsavename, 'LR')
    ROIs = {'Left', 'Right'};
    ROInames = ROIs;
elseif matches(ROIsavename, 'BigROI')
    ROIs = {'VisualROI_R','SensoryROI_R','MotorROI_R','RetrosplenialROI_R',...
        'VisualROI_L','SensoryROI_L','MotorROI_L','RetrosplenialROI_L'};
    ROInames = {'Visual Right', 'Sensory Right', 'Motor Right', 'Retrosplenial Right',...
        'Visual Left', 'Sensory Left', 'Motor Left', 'Retrosplenial Left'};
    % elseif matches(ROIsavename, 'WholeBrain')
    %     ROIs = {'WholeBrain'};
end

if ~exist('SaveDir', 'var')
    SaveDir =  '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist([SaveDir '/NVC/Sham_RS/'], 'dir')
    mkdir([SaveDir '/NVC/Sham_RS/'])
end

if ~exist('Acq', 'var')
    Acq = 'A3';
end

xaxeslimits = [-5 5];

%% Get specs and means/sems per acquisition
[allSpecs] = TableAllSpecs(RecordingOverview, Acq, ROIs, datanamefluocurves, NVCname);
allSpecs = sortrows(allSpecs, {'Group', 'Sex'});

if matches(ROIsavename, 'BigROI')
    for ind = 1:size(allSpecs,1)
        allSpecs.ROI{ind} = [allSpecs.ROI{ind}(1), allSpecs.ROI{ind}(end)];
    end
    save([SaveDir filesep 'NVC/Sham_RS/allSpecs.mat'], 'allSpecs'); %for stats later
end

[fluogroup, hbrgroup, hbogroup, means,sems] = ArrayAllCurves(RecordingOverview, Acq, ROIs, datanamefluocurves, NVCname);
[fluogroup_WB, hbrgroup_WB, hbogroup_WB, means_WB,sems_WB] = ArrayAllCurves(RecordingOverview, Acq, {'WholeBrain'}, datanamefluocurves, NVCname);


%% start plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Curves whole brain
f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout('flow');
x = linspace(-5, 10, 226);

%plot all mouselines
nexttile
title(['Per Mouse (N = ' num2str(size(RecordingOverview, 1)) ')'])
yhbo = hbogroup_WB.WholeBrain;
yhbr = hbrgroup_WB.WholeBrain;
yfluo = fluogroup_WB.WholeBrain;
yfluo = yfluo-1;

%right side
yyaxis right
plot(x, yhbo', '-', 'Color', [1 0 0 0.3])
hold on
plot(x, yhbr', '-', 'Color', [0 0 1 0.3]);
ylim([-2 2]);
% h_label = ylabel('\Delta \muM', 'interpreter', 'Tex', 'Rotation', 270);
h_label = ylabel('Hemodynamics (DmM)','Rotation', 270); % to work with illustrator (does not work with special characters)

ax = gca;
ax.YColor = 'red';
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
xlim(xaxeslimits);

%left side
yyaxis left
plot(x,yfluo, '-', 'Color', [0.4660 0.6740 0.1880 0.5]);

% f_label = ylabel('\Delta F/F');
f_label = ylabel('GCaMP (DF/F)'); % to work with illustrator
ylim([-0.05 0.05]); % centered at 0
ax = gca;
ax.YColor = [0.4660 0.6740 0.1880];
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
xlabel('Time (sec.)')
xlim([-5 10]);

leg1 = legend({'GCaMP', 'HbO', 'HbR'}, 'Orientation', 'Horizontal');
leg1.Location = 'southoutside';

% Plot average
nexttile
title('Average')
yhbo = means_WB.hbo.WholeBrain;
SEMhbo = sems_WB.hbo.WholeBrain;
% CIhbo = 1.96*SEMhbo/sqrt(size(RecordingOverview, 1));
yhbr = means_WB.hbr.WholeBrain;
SEMhbr = sems_WB.hbr.WholeBrain;
y = means_WB.fluo.WholeBrain-1;
SEM = sems_WB.fluo.WholeBrain;

%Right side
yyaxis right
plot(x, yhbo, 'Color', 'red', 'LineWidth', 2)
patch([x, fliplr(x)], [yhbo + SEMhbo fliplr(yhbo - SEMhbo)], 'r' ,'EdgeColor','none', 'FaceAlpha',0.25)
% patch([x, fliplr(x)], [yhbo + CIhbo fliplr(yhbo - CIhbo)], 'r' ,'EdgeColor','none', 'FaceAlpha',0.25)
hold on

% HbR
plot(x, yhbr, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2)
patch([x, fliplr(x)], [yhbr + SEMhbr fliplr(yhbr - SEMhbr)], 'b' ,'EdgeColor','none', 'FaceAlpha',0.25)

ylim([-2 2]);
% h_label = ylabel('\Delta \muM', 'interpreter', 'Tex', 'Rotation', 270);
h_label = ylabel('Hemodynamics (DmM)','Rotation', 270); % to work with illustrator

ax = gca;
ax.YColor = 'red';
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
xlim([-5 10]);

%Left side
yyaxis left
plot(x,y, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
patch([x, fliplr(x)], [y + SEM fliplr(y - SEM)], 'g' ,'EdgeColor','none', 'FaceAlpha',0.25)

% f_label = ylabel('\Delta F/F');
f_label = ylabel('GCaMP (DF/F)'); % to work with illustrator
ylim([-0.05 0.05]); % centered at 0
ax = gca;
ax.YColor = [0.4660 0.6740 0.1880];
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
xlabel('Time (sec.)')
xlim([-5 10]);

leg2 = legend({'GCaMP', 'SEM', 'HbO','SEM', 'HbR','SEM'}, 'Orientation', 'Horizontal');
leg2.Location = 'southoutside';
f.Position = [10 10 1800 800];
title(t, 'Whole Brain NVC')

%save
pause(0.5)
% saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_WholeBrain_Curves.tiff'], 'tiff');
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_WholeBrain_Curves.svg'], 'svg');
close(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% curves (per roi) - per mouse
f = figure('InvertHardcopy','off','Color',[1 1 1]);
if size(ROIs,2)>1
    t = tiledlayout(size(ROIs,2)/2,2);
else
    t = tiledlayout(1,1);
end
x = linspace(-5, 10, 226);

tileind = 2;
for roi = 1:size(ROIs, 2)
    nexttile(tileind)

    yhbo = hbogroup.(ROIs{roi});
    yhbr = hbrgroup.(ROIs{roi});
    yfluo = fluogroup.(ROIs{roi});
    yfluo = yfluo-1;

    yyaxis right
    plot(x, yhbo', '-', 'Color', [1 0 0 0.3])
    hold on
    plot(x, yhbr', '-', 'Color', [0 0 1 0.3]); %'blue')
    ylim([-2 2]);
    % h_label = ylabel('\Delta \muM', 'interpreter', 'Tex', 'Rotation', 270);

    ax = gca;
    ax.YColor = 'red';
    ax.XColor = 'k';
    set(ax, 'FontSize', 15, 'LineWidth', 2)
    xlim(xaxeslimits);

    %Left side
    % Fluo
    yyaxis left
    plot(x,yfluo, '-', 'Color', [0.4660 0.6740 0.1880 0.5]);

    % f_label = ylabel('\Delta F/F');
    ylim([-0.05 0.05]); % centered at 0
    ax = gca;
    ax.YColor = [0.4660 0.6740 0.1880];
    ax.XColor = 'k';
    set(ax, 'FontSize', 15, 'LineWidth', 2)
    xlim(xaxeslimits);

    xlabel('Time (sec)')
    title(ROInames{roi})

    if tileind == size(ROIs, 2)
        tileind = 1;
    else
        tileind = tileind+2;
    end
end

leg1 = legend({'GCaMP', 'HbO', 'HbR'}, 'Orientation', 'Horizontal');
leg1.Location = 'southoutside';
f.Position = [10 10 900 1000]; % if 2 wide
% f.Position = [10 10 1800 500]; % if 4 wide

% save
pause(0.5)
% saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_' ROIsavename '_IndividualCurves.tiff'], 'tiff');
% saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_' ROIsavename '_IndividualCurves.eps'], 'epsc');
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_' ROIsavename '_IndividualCurves.svg'], 'svg');
close(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% curves (per roi) - Averages
f = figure('InvertHardcopy','off','Color',[1 1 1]);
% t = tiledlayout(size(ROIs,2)/2, 2);
t = tiledlayout(size(ROIs,2)/4, 4);
x = linspace(-5, 10, 226);

tileind = 2;
for roi = 1:size(ROIs, 2)
    nexttile(tileind)
    title(ROInames{roi})

    yhbo = means.hbo.(ROIs{roi});
    SEMhbo = sems.hbo.(ROIs{roi});
    yhbr = means.hbr.(ROIs{roi});
    SEMhbr = sems.hbr.(ROIs{roi});
    y = means.fluo.(ROIs{roi});
    SEM = sems.fluo.(ROIs{roi});
    y = y-1; % to get centered around 0

    %% Modification : Clip data to x-axis limits
    idx = x >= xaxeslimits(1) & x <= xaxeslimits(2);
    x_clip = x(idx);
    y_clip = y(idx);
    yhbo_clip = yhbo(idx);
    yhbr_clip = yhbr(idx);
    SEM_clip = SEM(idx);
    SEMhbo_clip = SEMhbo(idx);
    SEMhbr_clip = SEMhbr(idx);

    % HbO
    yyaxis right
    plot(x_clip, yhbo_clip, 'Color', 'red', 'LineWidth', 2)
    patch([x_clip, fliplr(x_clip)], [yhbo_clip + SEMhbo_clip fliplr(yhbo_clip - SEMhbo_clip)], 'r' ,'EdgeColor','none', 'FaceAlpha',0.25)
    hold on
    % HbR
    plot(x_clip, yhbr_clip, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2)
    patch([x_clip, fliplr(x_clip)], [yhbr_clip + SEMhbr_clip fliplr(yhbr_clip - SEMhbr_clip)], 'b' ,'EdgeColor','none', 'FaceAlpha',0.25)
    ylim([-2 2]);
    % h_label = ylabel('\Delta \muM', 'interpreter', 'Tex', 'Rotation', 270);
    ax = gca;
    ax.YColor = 'red';
    ax.XColor = 'k';
    set(ax, 'FontSize', 15, 'LineWidth', 2)
        % set(ax,'yticklabel',[])
    xlim(xaxeslimits);
    % Fluo
    yyaxis left
    plot(x_clip, y_clip, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);
    patch([x_clip, fliplr(x_clip)], [y_clip + SEM_clip fliplr(y_clip - SEM_clip)], 'g' ,'EdgeColor','none', 'FaceAlpha',0.25)

    % f_label = ylabel('\Delta F/F');
    ylim([-0.05 0.05]); % centered at 0
    ax = gca;
    ax.YColor = [0.4660 0.6740 0.1880];
    ax.XColor = 'k';
    set(ax, 'FontSize', 15, 'LineWidth', 2)
    % set(ax,'yticklabel',[])
    xlabel('Time (sec)')
    xlim(xaxeslimits);

    if tileind == size(ROIs, 2)
        tileind = 1;
    else
        tileind = tileind+2;
    end
end

leg1 = legend({'GCaMP', 'SEM', 'HbO','SEM', 'HbR','SEM'}, 'Orientation', 'Horizontal');
leg1.Location = 'southoutside';
% f.Position = [10 10 1800 1000];
% f.Position = [10 10 900 1000];
f.Position = [10 10 1800 500];

% save
pause(0.5)
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_' ROIsavename '_AvCurves.svg'], 'svg');
close(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% boxplots
specs = {'GCaMPPeak', 'HbOPeak','HbRPeak', 'GCaMP_Increase', 'HbO_Increase', 'HbR_Decrease', ...
    'DelaySec', 'ResponseStrength', 'Strength_Increase'};
spectitles = {'GCaMP peak', 'HbO peak', 'HbR peak', 'GCaMP Increase', 'HbO Increase', 'HbR Decrease', 'Delay', 'Strength', 'Strength of Increase'};
% specylabels = {'\Delta F/F', '\Delta \muM', '\Delta \muM', '\Delta F/F Peak - Dip','\Delta \muM Peak - Dip','\Delta \muM Peak - Dip','Seconds', 'HbO/(GCaMP*100)', 'HbO Incr./(GCaMP Incr.*100)'};
specylabels = {'DF/F', 'DmM', 'DmM', 'DF/F Peak - Dip','DmM Peak - Dip','DmM Peak - Dip','Seconds', 'HbO/(GCaMP*100)', 'HbO Incr./(GCaMP Incr.*100)'};

f = figure('InvertHardcopy','off','Color',[1 1 1]);
% t = tiledlayout(size(specs,2),1);
t = tiledlayout('flow');

for indSpec = 1:size(specs,2)
    nexttile
    MakeBoxplotTile(allSpecs, 'ROI', specs{indSpec}, specylabels{indSpec});
    title(spectitles{indSpec})
end

f.Position = [10 10 1800 1000];
% title(t, ['Acquisition ' Acq])

% save
pause(0.5)
saveas(gcf, [SaveDir '/NVC/Sham_RS/' Acq '_' type '_' ROIsavename '_Boxplots.svg'], 'svg');
close(f)

end


function [allSpecs] = TableAllSpecs(RecordingOverview, Acquisition, ROInames, datanamefluocurves, NVCname)
for ind = 1:size(RecordingOverview, 1) %go per mouse
    Mouse = RecordingOverview.Mouse{ind};

    eval(['DataFolder = [RecordingOverview.' Acquisition '{ind} filesep];']);
    DataFolder = DataFolder;
    SaveFolder = [RecordingOverview.SaveDirectory{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];

    if matches(Mouse, 'M23')
        continue %M23 has not enough fluorescence, so skip
    elseif matches(Mouse, 'M14') && matches(Acquisition, 'A3')
        continue %M14 has really weird outliers for hbo/hbr data
    elseif matches(Mouse, 'M32') && (matches(Acquisition, 'A2') || matches(Acquisition, 'A3'))
        continue %M32 has damaged window, hbo/hbr data looks very bizarre
    elseif ~exist([SaveFolder NVCname '.mat'], 'file')
    % else
        disp(Mouse)
        if matches(NVCname, 'NVC_ROI_unweighted')
            GetAverageCurves_unweighted(SaveFolder, datanamefluocurves); % if you didnt make the nvc curves yet
        else
            GetAverageCurves(SaveFolder, NVCname, datanamefluocurves); % if you didnt make the nvc curves yet
        end
        disp(['GetAverageCurves done for ' Mouse])
    end

    load([SaveFolder NVCname '.mat'], 'Specs')
    % patch for when dips and peaks were not in GetAverageCurves yet:
    if ~sum(contains(Specs.Properties.VariableNames, 'HbRDipAfter'))
        regions = Specs.ROI;
        Specs = [];
        load([SaveFolder NVCname '.mat'], 'fluocurves', 'hbocurves', 'hbrcurves');
        imfreq = 15;

        for indroi = 1:size(regions,1)
            [tempSpecs] = GetSpecsPatch(fluocurves, hbocurves, hbrcurves, ...
                regions{indroi}, 5*imfreq, imfreq);
            Specs = [Specs; tempSpecs];
        end
        save([SaveFolder NVCname '.mat'], 'Specs', "-append")
    end

    Specsmouse = Specs;

    for ind2 = 1:size(ROInames, 2)
        Specs = Specsmouse(matches(Specsmouse.ROI,ROInames{ind2}),:);
        Specs.Properties.VariableNames{1} = 'Mouse';
        Specs.Mouse = {Mouse};
        Specs.Group = RecordingOverview.Group(ind);
        Specs.Sex = RecordingOverview.Sex(ind);
        Specs.Combi = RecordingOverview.Combi(ind);
        Specs.ROI = cellstr(ROInames{ind2});

        if ~exist('allSpecs', 'var')
            allSpecs = Specs;
        else
            allSpecs = [allSpecs; Specs];
        end
    end
end
end

function [Specs] = GetSpecsPatch(fluo, hbo, hbr, roiname, timebefore, imfreq)
VarNames = {'ROI', 'DelaySec', 'DelayFrames', 'ResponseStrength', 'ResponseStrengthRelative'...
    'GCaMPDipBefore', 'GCaMPPeak', 'GCaMPDipAfter',...
    'HbODipBefore', 'HbOPeak', 'HbODipAfter', ...
    'HbRDipBefore', 'HbRPeak', 'HbRDipAfter'};
VarTypes =     {'cell', 'single', 'single', 'single', 'single', ...
    'single', 'single', 'single',...
    'single', 'single', 'single',...
    'single', 'single', 'single'};
Specs = table('Size', size(VarNames), 'VariableNames', VarNames, ...
    'VariableTypes', VarTypes);


Specs.ROI = cellstr(roiname);

eval(['[maxfluo, indfluo] = findpeaks(fluo.' roiname '(timebefore:end));']) %find first peak after detected activation
if isempty(maxfluo) %if there's only nan
    Specs.DelaySec = NaN;
    Specs.DelayFrames = NaN;
    Specs.ResponseStrength = NaN;
    return
end

% GCaMPPeak
indfluo = timebefore + indfluo(1) - 1;
maxfluo = maxfluo(1);
Specs.GCaMPPeak = maxfluo;
% GCaMPDipBefore & GCaMPDipAfter
indDips = islocalmin(fluo.(roiname));
indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
Specs.GCaMPDipBefore = fluo.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
Specs.GCaMPDipAfter = fluo.(roiname)(indDipAfter);

%HbOPeak
eval(['[maxhbo, indhbo] = findpeaks(hbo.' roiname '(indfluo:end));'])
indhbo = indfluo + indhbo(1) - 1;
maxhbo = maxhbo(1);
Specs.HbOPeak = maxhbo;
% HbODipBefore & HbODipAfter
indDips = islocalmin(hbo.(roiname));
indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
Specs.HbODipBefore = hbo.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
Specs.HbODipAfter = hbo.(roiname)(indDipAfter);

%HbRPeak
[maxhbr, ~] = findpeaks(hbr.(roiname)(indfluo:end));
maxhbr = maxhbr(1);
Specs.HbRPeak = maxhbr;
% HbRDipBefore & HbRDipAfter
indDips = islocalmin(hbr.(roiname));
indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
Specs.HbRDipBefore = hbr.(roiname)(indDipBefore);
indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
Specs.HbRDipAfter = hbr.(roiname)(indDipAfter);

%Delays
Specs.DelaySec = (indhbo - indfluo)/imfreq; %in seconds
Specs.DelayFrames = indhbo - indfluo;
maxfluo = (maxfluo-1) *100; % in percentage

% Response Strengths
Specs.ResponseStrength = maxhbo / maxfluo;
increasefluo = Specs.GCaMPPeak - Specs.GCaMPDipBefore;
increasehbo = Specs.HbOPeak - Specs.HbODipBefore;
Specs.ResponseStrengthRelative = increasehbo/increasefluo;

end

function [fluogroup, hbrgroup, hbogroup, Means,SEMs] = ArrayAllCurves(RecordingOverview, Acquisition, ROInames, datanamefluocurves, NVCname)

% Get right size for curves by loading the first one
eval(['DataFolder = [RecordingOverview.' Acquisition '{1} filesep];']);
SaveFolder = [RecordingOverview.SaveDirectory{1} filesep ...
    RecordingOverview.Mouse{1} filesep DataFolder(end-5:end) 'CtxImg' filesep];
load([SaveFolder NVCname '.mat'], 'fluocurves')

curvelength = size(fluocurves.WholeBrain, 2);
clear fluocurves SaveFolder DataFolder

OverviewGroup = RecordingOverview;
fluogroup = [];
hbogroup = [];
hbrgroup = [];

for indMouse = 1:size(OverviewGroup, 1) %go per mouse
    Mouse = OverviewGroup.Mouse{indMouse};

    eval(['DataFolder = [OverviewGroup.' Acquisition '{indMouse} filesep];']);
    SaveFolder = [OverviewGroup.SaveDirectory{indMouse} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];

    if matches(Mouse, 'M23')
        continue
    elseif ~exist([SaveFolder NVCname '.mat'], 'file')
        disp(Mouse)
        if matches(NVCname, 'NVC_ROI_unweighted')
            GetAverageCurves_unweighted(SaveFolder, datanamefluocurves); % if you didnt make the nvc curves yet
        else
            GetAverageCurves(SaveFolder, datanamefluocurves); % if you didnt make the nvc curves yet
        end
        disp(['GetAverageCurves done for ' Mouse])
    end

    load([SaveFolder NVCname '.mat'], 'fluocurves', 'hbocurves', 'hbrcurves')

    for indROI = 1:size(ROInames, 2)

        eval(['fluogroup.' ROInames{indROI} '(indMouse,:) = fluocurves.' ROInames{indROI} ';'])
        eval(['hbogroup.' ROInames{indROI} '(indMouse,:) = hbocurves.' ROInames{indROI} ';'])
        eval(['hbrgroup.' ROInames{indROI} '(indMouse,:) = hbrcurves.' ROInames{indROI} ';'])
    end
end


for indROI = 1:size(ROInames, 2)
    eval(['Means.fluo.' ROInames{indROI} '= mean(fluogroup.' ROInames{indROI} ', 1, ''omitnan'');'])
    eval(['SEMs.fluo.' ROInames{indROI} ' = std(fluogroup.' ROInames{indROI} ', 0, 1,''omitnan'')/sqrt(size(fluogroup.' ROInames{indROI} ', 1));'])

    eval(['Means.hbo.' ROInames{indROI} '= mean(hbogroup.' ROInames{indROI} ', 1, ''omitnan'');'])
    eval(['SEMs.hbo.' ROInames{indROI} ' = std(hbogroup.' ROInames{indROI} ', 0, 1,''omitnan'')/sqrt(size(hbogroup.' ROInames{indROI} ', 1));'])

    eval(['Means.hbr.' ROInames{indROI} '= mean(hbrgroup.' ROInames{indROI} ', 1, ''omitnan'');'])
    eval(['SEMs.hbr.' ROInames{indROI} ' = std(hbrgroup.' ROInames{indROI} ', 0, 1,''omitnan'')/sqrt(size(hbrgroup.' ROInames{indROI} ', 1));'])
end

end


function MakeBoxplotTile(overviewtable, x, y, y_label)

if ~exist('x','var')
    x = 'ROI';
end

%name the columns in your table that you want to use x and y
temp = find(matches(overviewtable.Properties.VariableNames, x));
overviewtable.Properties.VariableNames{temp} = 'x';
overviewtable.x = categorical(overviewtable.x);
overviewtable.idx = grp2idx(overviewtable.x); %this gives the groups based on alphabet, so sort the ROI labels as well:
labels = cellstr(unique(overviewtable.x))';

if ~exist('y', 'var') && sum(contains(overviewtable.Properties.VariableNames, 'y')) == 0
    disp('Cannot make boxplot tile without Y variable')
    return
end

% Get the differences in peak/dip
if matches(y, 'GCaMP_Increase') %increase for gcamp/hbo, decrease for hbr
    overviewtable.y = overviewtable.GCaMPPeak - overviewtable.GCaMPDipBefore;
elseif matches(y, 'HbO_Increase')
    overviewtable.y = overviewtable.HbOPeak - overviewtable.HbODipBefore;
elseif matches(y, 'HbR_Decrease')
    overviewtable.y = overviewtable.HbRPeak - overviewtable.HbRDipAfter;
elseif matches(y, 'GCaMPPeak')
    overviewtable.y = overviewtable.GCaMPPeak - 1;
elseif matches(y, 'Strength_Increase')
    overviewtable.y = (overviewtable.HbOPeak - overviewtable.HbODipBefore)./...
        ((overviewtable.GCaMPPeak - overviewtable.GCaMPDipBefore)*100);
else
    temp = find(matches(overviewtable.Properties.VariableNames, y));
    overviewtable.Properties.VariableNames{temp} = 'y';
end


% BOXPLOT
%dont display outliers because we will do scatter that will show them
b = boxchart(overviewtable.idx, overviewtable.y, 'LineWidth', 2, 'MarkerStyle', 'none');
hold on
% hs = scatter(overviewtable.idx, overviewtable.y, 70, 'filled', 'jitter','on','JitterAmount',0.05);
hs = scatter(overviewtable.idx, overviewtable.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);

% MAKE PRETTY
% hs.MarkerFaceColor = [0 0 0];
hs.MarkerFaceColor = [0 0.4470 0.7410]; %blue because it's sham
hs.MarkerFaceAlpha = 0.5;
xticks(1:length(labels));
xticklabels(labels);
xlim([0.2 length(labels)+0.7])

axes1 = b.Parent;
hold(axes1,'on');
set(axes1, 'FontSize', 15, 'LineWidth', 2)
ylabel(y_label)

end

