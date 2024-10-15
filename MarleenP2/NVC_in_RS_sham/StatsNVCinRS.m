% load('/media/mbakker/GDrive/P2/GCaMP/NVC/Sham_RS/allSpecs.mat', 'allSpecs');


function StatsNVCinRS(SaveDir, allSpecs, type, Acq)

if ~exist('Acq', 'var')
    Acq = 'A1';
end

if ~exist('SaveDir', 'var')
    SaveDir = '/media/mbakker/GDrive/P2/GCaMP/';
elseif ~matches(SaveDir(end), filesep)
    SaveDir = [SaveDir filesep];
end

if ~exist('allSpecs', 'var')
    load([SaveDir 'NVC/Sham_RS/allSpecs_' type '_' Acq '.mat'], 'allSpecs');
end

specs = {'GCaMPPeak', 'HbOPeak','HbRPeak', 'GCaMP_Increase', 'HbO_Increase', 'HbR_Decrease', ...
    'DelaySec', 'ResponseStrength', 'Strength_Increase'};
% spectitles = {'GCaMP peak', 'HbO peak', 'HbR peak', 'GCaMP Increase', 'HbO Increase', 'HbR Decrease', 'Delay', 'Strength', 'Strength of Increase'};

% ROIs = {'VisualROI_R','SensoryROI_R','MotorROI_R','RetrosplenialROI_R',...
%         'VisualROI_L','SensoryROI_L','MotorROI_L','RetrosplenialROI_L'};
% ROIs = {'VR','SR','MR','RR','VL','SL','ML','RL'};

%% make good table
% NEW
% allSpecs.ResponseStrength = allSpecs.HbOPeak./allSpecs.GCaMPPeak;

% add extra columns
allSpecs.GCaMP_Increase = allSpecs.GCaMPPeak - allSpecs.GCaMPDipBefore;
allSpecs.HbO_Increase = allSpecs.HbOPeak - allSpecs.HbODipBefore;
allSpecs.HbR_Decrease = allSpecs.HbRPeak - allSpecs.HbRDipAfter;
allSpecs.Strength_Increase = allSpecs.HbO_Increase./(allSpecs.GCaMP_Increase*100);
% NEW
% allSpecs.Strength_Increase = allSpecs.HbO_Increase./(allSpecs.GCaMP_Increase);
allSpecs.GCaMPPeak = allSpecs.GCaMPPeak - 1;

% Move things around and clean up
allSpecs = movevars(allSpecs, "GCaMPPeak", "Before", "DelaySec");
allSpecs = movevars(allSpecs, "HbOPeak", "Before", "DelaySec");
allSpecs = movevars(allSpecs, "HbRPeak", "Before", "DelaySec");
allSpecs = movevars(allSpecs, "GCaMP_Increase", "Before", "DelaySec");
allSpecs = movevars(allSpecs, "HbO_Increase", "Before", "DelaySec");
allSpecs = movevars(allSpecs, "HbR_Decrease", "Before", "DelaySec");
allSpecs = removevars(allSpecs, "DelayFrames");
allSpecs = movevars(allSpecs, "Strength_Increase", "Before", "ResponseStrengthRelative");
allSpecs = removevars(allSpecs, ["ResponseStrengthRelative","GCaMPDipBefore","GCaMPDipAfter","HbODipBefore","HbODipAfter","HbRDipBefore","HbRDipAfter","Group","Sex","Combi"]);

% Make table per spec
nrofmice = size(unique(allSpecs.Mouse),1);
allSpecs = rows2vars(allSpecs);
temp = repmat(allSpecs(:,1), nrofmice,1);
ANOVAtable_all = allSpecs(:,1:9);
warning off
ANOVAtable_all(1:size(temp, 1),1) = temp;
warning on
clear temp

for indmice = 1:nrofmice-1
    ANOVAtable_all(indmice*11+1:indmice*11+11, 2:9) = ...
        allSpecs(1:11, indmice*8+2:indmice*8+9);
end

%change variable names
for ind = 2:9
    ANOVAtable_all.Properties.VariableNames(ind) = table2cell(ANOVAtable_all(11,ind));
end
Mouse = ANOVAtable_all(matches(ANOVAtable_all.OriginalVariableNames, 'Mouse'),2);
Mouse.Properties.VariableNames(1) = "Mouse";

ANOVAtable_all = movevars(ANOVAtable_all, "ML", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "MR", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "RL", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "RR", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "SL", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "SR", "Before", "VR");
ANOVAtable_all = movevars(ANOVAtable_all, "VL", "Before", "VR");

ROIs = ANOVAtable_all.Properties.VariableNames(2:9);

% per spec:
ANOVAtable_spec = [];
for ind = 1:size(specs,2)
    temp = [ANOVAtable_all(matches(ANOVAtable_all.OriginalVariableNames, specs{ind}),:), Mouse];
    for indname = 1:size(ROIs,2)
        temp.(ROIs{indname}) = cell2mat(temp.(ROIs{indname}));
    end
    ANOVAtable_spec.(specs{ind}) = temp;
end
clear Mouse allSpecs ind indmice



%% start statistics
sz = [length(specs) 4];
varTypes = ["string", "string", "double", "double"];
varNames = ["Spec", "Testtype", "p", "q"];
statsoutcome = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

statsoutcome.Spec(:) = specs;
% go per spec:
meas = table([1 2 3 4 5 6 7 8]', 'VariableNames', {'ROI'});

for indspec = 1:size(specs,2)

    % check assumptions anova:
    % - continuous scale for data: true
    % Note: delay is continuous but on a short scale

    % disp(specs{indspec})
    spectable = ANOVAtable_spec.(specs{indspec});
    specarray = table2array(spectable(:,2:9));

    % - normal distribution of data
    warning off
    [h1, ~] = adtest(reshape(specarray, 1, []));
    warning onn
    % f1 = figure;
    % histogram(specarray);
    % f2 = figure;
    % qqplot(reshape(specarray, 1, []));
    % close(f1, f2)

    rm = fitrm(spectable, "ML-VR~1", 'WithinDesign',meas);

    % - sphericity
    maunchlytbl = mauchly(rm);
    maunchly = maunchlytbl.pValue;
    % ranovatbl = ranova(rm); %note: (intercept):ROI sign means some ROI in the group differ sign.
    ranovatbl = ranova(rm, "WithinModel", "ROI");

    if h1 == 0 && maunchly > 0.05
        % repeated measures anova, not corrected
        statsoutcome.p(indspec) = ranovatbl.pValue(3);
        statsoutcome.Testtype(indspec) = {'ranova'};

    elseif h1 == 0 && maunchly < 0.05
        % repeated measures anova, epsilon correction
        % take greenhouse geisser because it's more strict, lower bounds
        % has been critized a lot (seen on internet)
        statsoutcome.p(indspec) = ranovatbl.pValueGG(3);
        statsoutcome.Testtype(indspec) = {'ranova-GG'};

    else
        clear rm maunchly
        % Friedman test: assumptions are
        % - measured on 3 or more occasions (ROI)
        % - group is random sample of population
        % - dependent variable is continuous

        %friedman: better than kruskal wallis becuase it's for repeated
        %measures specifically

        [p, ~, ~] = friedman(specarray, 1, 'off');
        statsoutcome.p(indspec) = p;
        statsoutcome.Testtype(indspec) = {'friedman'};

    end

end

%% FDR
statsoutcome.q(:) = mafdr(statsoutcome.p,'BHFDR', 'true');
save([SaveDir 'NVC/Sham_RS/Statsoutcome_' type '_' Acq '.mat'], 'statsoutcome');

%% post hoc tests
% Save for article
sz = [(length(ROIs)^2-length(ROIs))/2 length(specs)+2];
varTypes = ["string" "string" repmat("double", 1, length(specs))];
p_table = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames',["ROI_1" "ROI_2" specs]);
cld_table = table;
cld_table.ROI = ROIs';

% give right names for roi_1 and roi_2
indstart = 1;
for ind = 1:length(ROIs)
    % roi_1
    indend = indstart+length(ROIs)-1-ind;
    p_table.ROI_1(indstart:indend) = repmat(ROIs(ind), length(ROIs)-ind, 1);

    %roi_2
    p_table.ROI_2(indstart:indend) = ROIs(ind+1:length(ROIs));

    indstart = indend+1;
end

for indspec = 1:size(specs,2)
    if statsoutcome.q(indspec) > 0.05
        % disp([specs{indspec} ' has no sign. difference between brain regions'])
        continue
    end

    % disp(specs{indspec})
    spectable = ANOVAtable_spec.(specs{indspec});
    specarray = table2array(spectable(:,2:9));

    switch statsoutcome.Testtype(indspec)
        case 'ranova'
            % ranovatbl = ranova(rm);
            % multcompare(rm, 'ROI');
            % dont have any normal ranova, make code if it turns out you do
            % disp('ranova without GG - check if went right, make code if yes')

            rm = fitrm(spectable, "ML-VR~1", 'WithinDesign',meas);
            ranovatbl = ranova(rm, 'WithinModel', 'ROI');
            [c56] = multcompare(rm, 'ROI');

            % remove double p values
            c = table;
            for ind = 1:size(c56,1)
                if sum(c56.ROI_1(1:ind) == c56.ROI_2(ind)) == 0
                    c = [c; c56(ind,:)];
                end
            end

            p_table.(specs{indspec}) = c.pValue;
            cld_temp = cld(c56);
            cld_table.(specs{indspec}) = cld_temp.letter;

        case  "ranova-GG"
            rm = fitrm(spectable, "ML-VR~1", 'WithinDesign',meas);
            ranovatbl = ranova(rm, 'WithinModel', 'ROI');
            [c56] = multcompare(rm, 'ROI');

            % remove double p values
            c = table;
            for ind = 1:size(c56,1)
                if sum(c56.ROI_1(1:ind) == c56.ROI_2(ind)) == 0
                    c = [c; c56(ind,:)];
                end
            end

            p_table.(specs{indspec}) = c.pValue;
            cld_temp = cld(c56);
            cld_table.(specs{indspec}) = cld_temp.letter;

        case 'friedman'
            [~, ~, stats] = friedman(specarray, 1, 'off');
            [c, ~, ~, ~] = multcompare(stats, 'Display', 'off'); %multcompare default is tukey-kramer

            c = array2table(c, "VariableNames",...
                ["ROI_1", "ROI_2", "Lower Limit", "ROI_1-ROI_2", "Upper Limit", "pValue"]);
            %note: m/means and difference between groups is based on rank

            p_table.(specs{indspec}) = c.pValue;
            c2 = c;
            c2.Properties.VariableNames(1) = "ROI_21";
            c2.Properties.VariableNames(2) = "ROI_1";
            c2.Properties.VariableNames(1) = "ROI_2";
            c2 = movevars(c2, 1, "Before", 3);
            c56 = [c; c2];
            c56 = sortrows(c56);
            cld_temp = cld(c56);
            cld_table.(specs{indspec}) = cld_temp.letter;

    end

    posthocoutcome.(specs{indspec}) = c;

end

%% save 
save([SaveDir 'NVC/Sham_RS/Posthoc_outcome_' type '_' Acq '.mat'], 'posthocoutcome', 'cld_table');
writetable(p_table, [SaveDir 'NVC/Sham_RS/Posthoc_outcome_' type '_' Acq '.xlsx'])
writetable(statsoutcome, [SaveDir 'NVC/Sham_RS/Stats_outcome_' type '_' Acq '.xlsx'])

end
