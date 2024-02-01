% Assumptions repeated measurement anova:
% 1. normal distribution -- Anderson-Darling
% 2. experimental group had 3 or more levels
% 3. no outliers
% 4. Sphericity -- Mauchly's test


function results = MovementStats(overviewtable)

% remove outlier(s)
% overviewtable.Movement(overviewtable.Movement>2000) = NaN;
outliers = isoutlier(overviewtable.Movement);
overviewtable.Movement(outliers) = NaN;
disp([num2str(sum(outliers)) ' outliers removed'])

% check normal distribution
Acquisitions = unique(overviewtable.Acquisition);
groups = unique(overviewtable.Combi);

for indacq = 1:size(Acquisitions, 2)
    groupdataAcq = overviewtable(overviewtable.Acquisition == Acquisitions(indacq),:);

    for indgroup = 1:size(groups, 1)
        groupdata = groupdataAcq(groupdataAcq.Combi == groups(indgroup),:);
        if adtest(groupdata.Movement) ~= 0
            disp(['Acq ' num2str(indacq) ' group ' char(groups(indgroup)) ' not normal distribution. Find different stats test.'])
            return
        end

    end
end

% get right table
Acq1 = overviewtable(overviewtable.Acquisition == 'A1',:);
Acq1.Properties.VariableNames(6) = "MovementA1";
Acq1 = removevars(Acq1, "Acquisition");

Acq2 = overviewtable(overviewtable.Acquisition == 'A2',{'Mouse' 'Movement'}');
Acq2.Properties.VariableNames(2) = "MovementA2";

Acq3 = overviewtable(overviewtable.Acquisition == 'A3',{'Mouse' 'Movement'});
Acq3.Properties.VariableNames(2) = "MovementA3";

rmtable = join(Acq1, Acq2);
rmtable = join(rmtable, Acq3);
clear Acq1 Acq2 Acq3 

% check for sphericity
Meas = dataset([1 2 3]', 'VarNames', {'Measurements'});
rm = fitrm(rmtable, "MovementA1-MovementA3~Combi", "WithinDesign",Meas);
mauchlyresults = mauchly(rm);

if mauchlyresults.pValue < 0.05
    disp('Assumption for sphericity not met, find different stats test or do epsilon correction.')
    return
end

% Do the actualy repeated measures anova
results = ranova(rm);

if results.pValue(1)<0.05
    disp(['Effect of week of testing on movement significant! p = ' num2str(results.pValue(1))])
end
if results.pValue(1)<0.05
    disp(['Interaction between week of testing and group significant! p = ' num2str(results.pValue(2))])
end

end

