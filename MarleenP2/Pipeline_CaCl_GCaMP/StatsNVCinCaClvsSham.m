% ROIs can be 'R' or 'BigROI'

function StatsNVCinCaClvsSham(ROIs, allSpecs)

if ~exist('ROIs', 'var')
    ROIs = 'R';
end

if ~exist("allSpecs", 'var')
    load(['/media/mbakker/GDrive/P2/GCaMP/NVC/allSpecs_' ROIs '.mat'], 'allSpecs');
end


specs = {'GCaMPPeak', 'HbOPeak','HbRPeak', 'GCaMP_Increase', 'HbO_Increase', 'HbR_Decrease', ...
    'DelaySec', 'ResponseStrength', 'Strength_Increase'};
Acqs = {'A1', 'A2', 'A3'};


%% make good tables
for indacq = 1:size(Acqs,2)
    Acq = Acqs{indacq};

    % add extra columns
    allSpecs.(Acq).GCaMP_Increase = allSpecs.(Acq).GCaMPPeak - allSpecs.(Acq).GCaMPDipBefore;
    allSpecs.(Acq).HbO_Increase = allSpecs.(Acq).HbOPeak - allSpecs.(Acq).HbODipBefore;
    allSpecs.(Acq).HbR_Decrease = allSpecs.(Acq).HbRPeak - allSpecs.(Acq).HbRDipAfter;
    allSpecs.(Acq).Strength_Increase = allSpecs.(Acq).HbO_Increase./(allSpecs.(Acq).GCaMP_Increase*100);
    allSpecs.(Acq).GCaMPPeak = allSpecs.(Acq).GCaMPPeak - 1;
end


% get table per spec
mice = unique(allSpecs.A1.Mouse);

for indspec = 1:size(specs, 2)
    T1 = table(allSpecs.A1.Mouse, allSpecs.A1.Group, allSpecs.A1.(specs{indspec}),...
        'VariableNames', ["Mice", "Group", "A1"]);

    T2 = table(allSpecs.A2.Mouse, allSpecs.A2.Group, allSpecs.A2.(specs{indspec}),...
        'VariableNames', ["Mice", "Group", "A2"]);
    missingmice = T1(~matches(allSpecs.A1.Mouse, T2.Mice),1:2);
    missingmice.A2 = NaN(size(missingmice,1),1);
    T2 = [T2; missingmice];

    T3 = table(allSpecs.A3.Mouse, allSpecs.A3.Group, allSpecs.A3.(specs{indspec}),...
        'VariableNames', ["Mice", "Group", "A3"]);
    missingmice = T1(~matches(allSpecs.A1.Mouse, T3.Mice),1:2);
    missingmice.A3 = NaN(size(missingmice,1),1);
    T3 = [T3; missingmice];

    T = join(T1,T2);
    T = join(T,T3);
    spectable.(specs{indspec}) = T;
end
clear T T1 T2 T3 missingmice indacq indspec

%% Stats 
sz = [length(specs) 4];
varTypes = ["string", "string", "double", "double"];
varNames = ["Spec", "Testtype", "p", "q"];
statsoutcome = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
statsoutcome.Spec(:) = specs;

% per spec
for indspec = 1:size(specs,2)

    spectable.(specs{indspec}).Group = string(spectable.(specs{indspec}).Group);
   
    % - normal distribution of data
    specarray = table2array(spectable.(specs{indspec})(:,3:5));
    [h1, ~] = adtest(reshape(specarray, 1, []));
    f1 = figure;
    histogram(specarray);
    close(f1)
    
    % fit model
    Time = [1 2 3]';
    rm = fitrm(spectable.(specs{indspec}), 'A1-A3 ~ Group', 'WithinDesign', Time);
    % ranovatbl = ranova(rm);
    ranovatbl = ranova(rm, "WithinModel", 'Time');

    % - sphericity
    maunchlytbl = mauchly(rm);
    maunchly = maunchlytbl.pValue;

   if h1 == 0 && maunchly > 0.05
        % repeated measures anova, not corrected
        statsoutcome.p(indspec) = ranovatbl.pValue(2);
        statsoutcome.Testtype(indspec) = {'ranova'};

    elseif h1 == 0 && maunchly < 0.05
        % repeated measures anova, epsilon correction
        statsoutcome.p(indspec) = ranovatbl.pValueGG(2);
        statsoutcome.Testtype(indspec) = {'ranova-GG'};

    else
        clear rm maunchly
        % cannot do Friedman because sample sizes are different (cacl 14
        % sham 17 or 15 if you discount the nans for a2 and a3)
        % Mack-Skillings test 
        % https://github.com/thomaspingel/mackskill-matlab
        % take hollander and wolfe example that is in code
        % Get table in right way:
        T1 = spectable.(specs{indspec})(:,1:3);
        T1.Properties.VariableNames(3) = "y";
        T1.Acq = repmat("A1", size(T1, 1), 1);

        T2 = [spectable.(specs{indspec})(:,1:2) spectable.(specs{indspec})(:,4)];
        T2.Properties.VariableNames(3) = "y";
        T2.Acq = repmat("A2", size(T2, 1), 1);

        T3 = [spectable.(specs{indspec})(:,1:2) spectable.(specs{indspec})(:,5)];
        T3.Properties.VariableNames(3) = "y";
        T3.Acq = repmat("A3", size(T3, 1), 1);      

        T = [T1; T2; T3];

        [mackskill_p, ~] = mackskill(T.y, T.Group, T.Acq);
        statsoutcome.p(indspec) = mackskill_p;
        statsoutcome.Testtype(indspec) = {'MackSkillings'};

   end
end

%% FDR
disp("A3")
statsoutcome.q(:) = mafdr(statsoutcome.p,'BHFDR', 'true');
save('/media/mbakker/GDrive/P2/GCaMP/NVC/CaClvsSham/Statsoutcome.mat', 'statsoutcome');

%% post hoc tests
% Did not find sign. difference between groups, so no post-hoc was done. 

end







