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
% - 'R', right cacl vs sham
% - 'BigROI', which takes all areas (visual, sensory etc.) except auditory
% - 'WholeBrain', no ROI, just the brain in general

function plotNVCAllAcq(ROIsavename, type, SaveDir)

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');

if ~exist('type','var') || matches(type, 'normal')
    type = 'normal';
    datanamefluocurves = 'hemoCorr_fluo';
    NVCname = 'NVC_ROI';
% elseif matches(type, 'OutliersRemoved')
%     datanamefluocurves = 'hemoCorr_fluo';
%     NVCname = 'NVC_ROI_OL_removed';
elseif matches(type, 'nofilt')
    datanamefluocurves = 'fluo_nofilt';
    NVCname = ['NVC_ROI_' datanamefluocurves];
elseif matches(type, 'unweighted')
    NVCname = 'NVC_ROI_unweighted';
end

if ~exist('ROIsavename', 'var')
    ROIsavename = 'LR';
end

if matches(ROIsavename, 'LR')
    ROIs = {'Left', 'Right'};
    ROInames = ROIs;
elseif matches(ROIsavename, 'R')
    ROIs = {'Right'};
    ROInames = ROIs;
elseif matches(ROIsavename, 'BigROI')
    ROIs = {'VisualROI_R','SensoryROI_R','MotorROI_R','RetrosplenialROI_R',...
        'VisualROI_L','SensoryROI_L','MotorROI_L','RetrosplenialROI_L'};
    ROInames = {'Visual Right', 'Sensory Right', 'Motor Right', 'Retrosplenial Right',...
        'Visual Left', 'Sensory Left', 'Motor Left', 'Retrosplenial Left'};
elseif matches(ROIsavename, 'WholeBrain')
    ROIs = {'WholeBrain'};
    ROInames = ROIs;
end

if ~exist('SaveDir', 'var')
    SaveDir =  '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist([SaveDir '/NVC/Specs/'], 'dir')
    mkdir([SaveDir '/NVC/Specs/'])
end

xaxeslimits = [-5 5];

Acquisitions = {'A1', 'A2', 'A3'};
% Acquisitions = {'A1'};

specs = {'GCaMPPeak', 'HbOPeak','HbRPeak', 'GCaMP_Increase', 'HbO_Increase', 'HbR_Decrease', ...
    'DelaySec', 'ResponseStrength', 'Strength_Increase'};
spectitles = {'GCaMP peak', 'HbO peak', 'HbR peak', 'GCaMP Increase', 'HbO Increase', 'HbR Decrease', 'Delay', 'Strength', 'Strength of Increase'};
specylabels = {'DF/F', 'DmM', 'DmM', 'DF/F Peak - Dip','DmM Peak - Dip','DmM Peak - Dip','Seconds', 'HbO/(GCaMP*100)', 'HbO Incr./(GCaMP Incr.*100)'};


%% Choose grouping
Grouping = {'Group'};
[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

%% Get specs and means/sems per acquisition
allSpecs = [];
for indAcq = 1:size(Acquisitions, 2)
    Acq = Acquisitions{indAcq};

    [AcqSpecs] = TableAllSpecs(RecordingOverview, Acq, ROIs, datanamefluocurves, NVCname);
    AcqSpecs = sortrows(AcqSpecs, {'Group', 'Sex'});

    if matches(ROIsavename, 'BigROI')
        for ind = 1:size(AcqSpecs,1)
            AcqSpecs.ROI{ind} = [AcqSpecs.ROI{ind}(1), AcqSpecs.ROI{ind}(end)];
        end
    end
    allSpecs.(Acq) = AcqSpecs;

    % [Acqmeans,Acqsems] = ArrayAllCurves(Grouping, groups, RecordingOverview, Acq, ROIs, datanamefluocurves, NVCname);
    % means.(Acq) = Acqmeans;
    % sems.(Acq) = Acqsems;

    clear AcqSpecs Acqmeans Acqsems ind
end
clear Acq indAcq

if matches(ROIsavename, 'BigROI')
    save([SaveDir filesep 'NVC/allSpecs_BigROI.mat'], 'allSpecs'); %for stats later
elseif matches(ROIsavename, 'R')
    save([SaveDir filesep 'NVC/allSpecs_R.mat'], 'allSpecs'); %for stats later
end

%% Get ylim values
% make sure theyre the same for all 3 acquisitions so you can compare them
% to each other but still see all the individual points
ylims = zeros(size(specs,2),2);

for indspec = 1:size(specs, 2)
    temp = [];

    for indAcq = 1:size(Acquisitions, 2)
        Acq = Acquisitions{indAcq};

        if matches(specs{indspec}, 'GCaMP_Increase')
            allSpecs.(Acq).(specs{indspec}) = allSpecs.(Acq).GCaMPPeak - allSpecs.(Acq).GCaMPDipBefore;
        elseif matches(specs{indspec}, 'HbO_Increase')
            allSpecs.(Acq).(specs{indspec}) = allSpecs.(Acq).HbOPeak - allSpecs.(Acq).HbODipBefore;
        elseif matches(specs{indspec}, 'HbR_Decrease')
            allSpecs.(Acq).(specs{indspec}) = allSpecs.(Acq).HbRPeak - allSpecs.(Acq).HbRDipAfter;
        elseif matches(specs{indspec}, 'Strength_Increase')
            allSpecs.(Acq).(specs{indspec}) = allSpecs.(Acq).HbO_Increase ./...
                (allSpecs.(Acq).GCaMP_Increase .*100);
        end
        temp = [temp; allSpecs.(Acq).(specs{indspec})];
    end

    if matches(specs{indspec}, 'GCaMPPeak')
        temp = temp-1;
    end

    ylims(indspec,1) = min(temp, [], "all");
    ylims(indspec,2) = max(temp, [], "all");
end

%% start plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% boxplots
if length(ROIs) == 1
    f = figure('InvertHardcopy','off','Color',[1 1 1]);
    t = tiledlayout(3,3);

    specs3acqs = [];
    for indacq = 1:size(Acquisitions, 2)
        specsacq = allSpecs.(Acquisitions{indacq});
        temp = cell(size(specsacq, 1), 1);
        temp(:) = {Acquisitions{indacq}};
        specsacq.Acq = temp;
        specs3acqs = [specs3acqs; specsacq];
    end
    clear specsacq indacq temp

    for indSpec = 1:size(specs,2)
        nexttile
        MakeBoxplotTile(Grouping, groups, specs3acqs, 'Acq', specs{indSpec}, specylabels{indSpec}, ylims(indSpec,:))
        title(spectitles{indSpec})
    end
    f.Position = [10 10 1800 1000];
    title(t,'Right Hemisphere'); %change if you alter roi

    % save
    saveas(gcf, [SaveDir '/NVC/Specs/AllAcq_' type '_' ROIsavename '_Boxplots.svg'], 'svg');
    close(f)
else

    for indAcq = 1:size(Acquisitions,2)
        f = figure('InvertHardcopy','off','Color',[1 1 1]);
        % t = tiledlayout(size(Acquisitions,2), size(specs,2));
        t = tiledlayout(3,3);
        Acq = Acquisitions{indAcq};
        AcqSpecs = allSpecs.(Acq);

        for indSpec = 1:size(specs,2)
            nexttile
            MakeBoxplotTile(Grouping, groups, AcqSpecs, 'ROI', specs{indSpec}, specylabels{indSpec}, ylims(indSpec,:))
            title(spectitles{indSpec})
        end

        f.Position = [10 10 1800 1000];
        title(t,['Acquisition ' Acq]);

        % save
        saveas(gcf, [SaveDir '/NVC/Specs/' Acq '_' type '_' ROIsavename '_Boxplots.svg'], 'svg');
        close(f)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allSpecs] = TableAllSpecs(RecordingOverview, Acquisition, ROInames, datanamefluocurves, NVCname)

for ind = 1:size(RecordingOverview, 1) %go per mouse
    Mouse = RecordingOverview.Mouse{ind};

    % eval(['DataFolder = [RecordingOverview.' Acquisition '{ind} filesep];']);
    DataFolder = [RecordingOverview.(Acquisition){ind} filesep];
    % DataFolder = DataFolder;
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
        end
        disp(['GetAverageCurves done for ' Mouse])
    end

    load([SaveFolder NVCname '.mat'], 'Specs')

    % %temp
    % patch for when dips and peaks were not in GetAverageCurves yet:
    % if ~sum(contains(Specs.Properties.VariableNames, 'HbRDipAfter'))
    %     regions = Specs.ROI;
    %     Specs = [];
    %     load([SaveFolder NVCname '.mat'], 'fluocurves', 'hbocurves', 'hbrcurves');
    %     imfreq = 15;
    % 
    %     for indroi = 1:size(regions,1)
    %         [tempSpecs] = GetSpecsPatch(fluocurves, hbocurves, hbrcurves, ...
    %             regions{indroi}, 5*imfreq, imfreq);
    %         Specs = [Specs; tempSpecs];
    %     end
    %     save([SaveFolder NVCname '.mat'], 'Specs', "-append")
    % end

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


% function [Specs] = GetSpecsPatch(fluo, hbo, hbr, roiname, timebefore, imfreq)
% VarNames = {'ROI', 'DelaySec', 'DelayFrames', 'ResponseStrength', 'ResponseStrengthRelative'...
%     'GCaMPDipBefore', 'GCaMPPeak', 'GCaMPDipAfter',...
%     'HbODipBefore', 'HbOPeak', 'HbODipAfter', ...
%     'HbRDipBefore', 'HbRPeak', 'HbRDipAfter'};
% VarTypes =     {'cell', 'single', 'single', 'single', 'single', ...
%     'single', 'single', 'single',...
%     'single', 'single', 'single',...
%     'single', 'single', 'single'};
% Specs = table('Size', size(VarNames), 'VariableNames', VarNames, ...
%     'VariableTypes', VarTypes);
% 
% 
% Specs.ROI = cellstr(roiname);
% 
% [maxfluo, indfluo] = findpeaks(fluo.(roiname)(timebefore:end)); %find first peak after detected activation
% if isempty(maxfluo) %if there's only nan
%     Specs.DelaySec = NaN;
%     Specs.DelayFrames = NaN;
%     Specs.ResponseStrength = NaN;
%     return
% end
% 
% % GCaMPPeak
% indfluo = timebefore + indfluo(1) - 1;
% maxfluo = maxfluo(1);
% Specs.GCaMPPeak = maxfluo;
% % GCaMPDipBefore & GCaMPDipAfter
% indDips = islocalmin(fluo.(roiname));
% indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
% Specs.GCaMPDipBefore = fluo.(roiname)(indDipBefore);
% indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
% Specs.GCaMPDipAfter = fluo.(roiname)(indDipAfter);
% 
% %HbOPeak
% eval(['[maxhbo, indhbo] = findpeaks(hbo.' roiname '(indfluo:end));'])
% indhbo = indfluo + indhbo(1) - 1;
% maxhbo = maxhbo(1);
% Specs.HbOPeak = maxhbo;
% % HbODipBefore & HbODipAfter
% indDips = islocalmin(hbo.(roiname));
% indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
% Specs.HbODipBefore = hbo.(roiname)(indDipBefore);
% indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
% Specs.HbODipAfter = hbo.(roiname)(indDipAfter);
% 
% %HbRPeak
% [maxhbr, ~] = findpeaks(hbr.(roiname)(indfluo:end));
% maxhbr = maxhbr(1);
% Specs.HbRPeak = maxhbr;
% % HbRDipBefore & HbRDipAfter
% indDips = islocalmin(hbr.(roiname));
% indDipBefore = find(indDips(1:indfluo(1)), 1, 'last');
% Specs.HbRDipBefore = hbr.(roiname)(indDipBefore);
% indDipAfter = indDipBefore + find(indDips(indDipBefore+1:end), 1, 'first');
% Specs.HbRDipAfter = hbr.(roiname)(indDipAfter);
% 
% %Delays
% Specs.DelaySec = (indhbo - indfluo)/imfreq; %in seconds
% Specs.DelayFrames = indhbo - indfluo;
% maxfluo = (maxfluo-1) *100; % in percentage
% 
% % Response Strengths
% Specs.ResponseStrength = maxhbo / maxfluo;
% increasefluo = Specs.GCaMPPeak - Specs.GCaMPDipBefore;
% increasehbo = Specs.HbOPeak - Specs.HbODipBefore;
% Specs.ResponseStrengthRelative = increasehbo/increasefluo;
% 
% end
% 


% only  need this when you want to plot curves:
function [Means,SEMs] = ArrayAllCurves(Grouping, groups, RecordingOverview, Acquisition, ROInames, datanamefluocurves, NVCname)
% if matches(datanamefluocurves, 'hemoCorr_fluo')
%     NVCname = 'NVC_ROI';
% else
%     NVCname = ['NVC_ROI_' datanamefluocurves];
% end

% Get right size for curves by loading the first one
eval(['DataFolder = [RecordingOverview.' Acquisition '{1} filesep];']);
SaveFolder = [RecordingOverview.SaveDirectory{1} filesep ...
    RecordingOverview.Mouse{1} filesep DataFolder(end-5:end) 'CtxImg' filesep];
load([SaveFolder NVCname '.mat'], 'fluocurves')

curvelength = size(fluocurves.WholeBrain, 2);
clear fluocurves SaveFolder DataFolder

% start building arrays
for indgroup = 1:size(groups,1) % Go per group
    group = groups{indgroup};
    eval(['OverviewGroup = RecordingOverview(RecordingOverview.' char(Grouping) ' == ''' group ''',:);'])
    fluogroup = NaN(size(OverviewGroup,1), curvelength, size(ROInames, 2));
    hbogroup = NaN(size(OverviewGroup,1), curvelength, size(ROInames, 2));
    hbrgroup = NaN(size(OverviewGroup,1), curvelength, size(ROInames, 2));


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
            eval(['fluogroup(indMouse,:,indROI) = fluocurves.' ROInames{indROI} ';'])
            eval(['hbogroup(indMouse,:, indROI) = hbocurves.' ROInames{indROI} ';'])
            eval(['hbrgroup(indMouse,:, indROI) = hbrcurves.' ROInames{indROI} ';'])
        end
    end

    % note to self, if you want a plot where all mice have a seperate line,
    % save the fluo/hbo/hbrgroup in the indROI forloop with the roinames
    % attached to them or something similar....
    for indROI = 1:size(ROInames, 2)
        x = fluogroup(:,:,indROI);
        eval(['Means.fluo.' group '.' ROInames{indROI} '= mean(x, 1, ''omitnan'');'])
        SEM = std(x, 0, 1,'omitnan')/sqrt(size(x, 1));
        eval(['SEMs.fluo.' group '.' ROInames{indROI} ' = SEM;'])

        x = hbogroup(:,:,indROI);
        eval(['Means.hbo.' group '.' ROInames{indROI} '= mean(x, 1, ''omitnan'');'])
        SEM = std(x, 0, 1,'omitnan')/sqrt(size(x, 1));
        eval(['SEMs.hbo.' group '.' ROInames{indROI} ' = SEM;'])

        x = hbrgroup(:,:,indROI);
        eval(['Means.hbr.' group '.' ROInames{indROI} '= mean(x, 1, ''omitnan'');'])
        SEM = std(x, 0, 1,'omitnan')/sqrt(size(x, 1));
        eval(['SEMs.hbr.' group '.' ROInames{indROI} ' = SEM;'])
    end

end
end

% Make sure in your overviewtable that you have:
% Mouse, some variable for x (usually ROI), Group, Sex, Combi (or at least the
% grouping that you want to use) and the y value that you want plotted.
%

function MakeBoxplotTile(Grouping, groups, overviewtable, x, y, y_label, ylimvalues)

if ~exist('x','var')
    x = 'ROI';
end

Grouping = char(Grouping);

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
eval(['b = boxchart(overviewtable.idx, overviewtable.y, ''GroupByColor'', overviewtable.' Grouping ', ''LineWidth'', 2, ''MarkerStyle'', ''none'');'])
hold on

% SCATTER
xaxisstep = 1/size(groups,1);
xaxisplacement = 1 - 0.5*xaxisstep - (size(groups,1)/2-1)*xaxisstep - xaxisstep;

for indroi = 1:length(labels)
    currentROI = overviewtable.x == labels{indroi}; %is called currentroi because it's usually grouped by ROI, but can be something else

    for indgroup = 1:size(groups,1)
        currentgroup = overviewtable.(Grouping) == groups{indgroup};
        % eval(['currentgroup = overviewtable.' Grouping ' == groups{indgroup};']);
        currentindex = currentgroup.*currentROI;
        overviewtable.idx(currentindex == 1) = xaxisplacement + xaxisstep*indgroup;
    end

    xaxisplacement = xaxisplacement + 1;
end
clear indroi indgroup currentgroup currentROI currentindex ind indmouse tablemouse xaxisplacement group idx Mousegroup toplot

% Scatter per group so you can vary colours
if matches(Grouping, 'Group') % hardcoded
    cacl = overviewtable(overviewtable.(Grouping)=='CaCl',:);
    hs = scatter(cacl.idx, cacl.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0.6350 0.0780 0.1840];
    hs.MarkerFaceAlpha = 0.5;
    clear hs cacl

    sham = overviewtable(overviewtable.(Grouping)=='Sham',:);
    hs = scatter(sham.idx, sham.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0 0.4470 0.7410];
    hs.MarkerFaceAlpha = 0.5;
    clear hs sham
else
    hs = scatter(overviewtable.idx, overviewtable.y, 20, 'filled', 'jitter','on','JitterAmount',0.15);
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.5;
end

xticks(1:length(labels));
xticklabels(labels);
xlim([0.2 length(labels)+0.7])

axes1 = b.Parent;
hold(axes1,'on');
set(axes1, 'FontSize', 15, 'LineWidth', 2)
% set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2);
% legend(b, 'Location', 'northeast', 'NumColumns', 2);
if exist('ylimvalues', 'var')
    ylim(ylimvalues);
end

b(1).SeriesIndex = 7;
b(2).SeriesIndex = 1;
if size(b,1)>2
    b(3).SeriesIndex = 2;
    b(4).SeriesIndex = 6;
end

ylabel(y_label)

end
