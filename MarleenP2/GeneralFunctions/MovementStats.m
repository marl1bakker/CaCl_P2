% Assumptions repeated measurement anova:
% 1. normal distribution -- Anderson-Darling
% 2. experimental group had 3 or more levels
% 3. no outliers
% 4. Sphericity -- Mauchly's test


function results = MovementStats(overviewtable)

% check normal distribution
warning('off','all');
normaldistr = adtest(overviewtable.Movement);
warning('on','all');

switch(normaldistr)
    %% normally distributed, do rmanova
    case 0
        
        % remove outlier(s)
        % overviewtable.Movement(overviewtable.Movement>2000) = NaN;
        outliers = isoutlier(overviewtable.Movement);
        overviewtable.Movement(outliers) = NaN;
        disp([num2str(sum(outliers)) ' outliers removed'])
        
        Acquisitions = unique(overviewtable.Acquisition);
        groups = unique(overviewtable.Combi);
        
        for indacq = 1:size(Acquisitions, 1)
            groupdataAcq = overviewtable(overviewtable.Acquisition == Acquisitions(indacq),:);
            
            for indgroup = 1:size(groups, 1)
                groupdata = groupdataAcq(groupdataAcq.Combi == groups(indgroup),:);
                if adtest(groupdata.Movement) ~= 0
                    disp(['Acq ' num2str(indacq) ' group ' char(groups(indgroup)) ' not normal distribution. Find different stats test.'])
                    results = [];
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
            results = [];
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
        
        %% not normally distributed
    case 1 
        %         %% not normally distributed, do friedmans
        %         %cant do friedman because sample sizes differ, do kruskall wallis
        
        %         % get good table
        %         groups = cellstr(unique(overviewtable.Combi));
        %         tablefriedman = table();
        %
        %         for indgroup = 1:size(groups, 1)
        %             group = overviewtable(overviewtable.Combi == groups{indgroup},:);
        %             group = sortrows(group, 'Mouse');
        %
        %             if size(group,1) < size(tablefriedman,1)
        %                 group.Movement(end+1:size(tablefriedman,1)) = NaN;
        %             elseif size(group,1) > size(tablefriedman,1) && indgroup ~= 1
        %                 rowdifference = size(group,1)-size(tablefriedman,1);
        %                 tablefriedman{end+1:end+rowdifference,:} = missing;
        %             end
        %
        %             eval(['tablefriedman.' groups{indgroup} ' = group.Movement;']);
        %         end
        %
        %         clear rowdifference group indgroup groups
        %
        %         % Friedmans test
        %         p = friedman(table2array(tablefriedman), 3); %reps is 3 for A1, A2, A3. Organisation of table is important. Did that above.
        
        %% kruskal wallis
        %go per acquisition
        Acquisitions = unique(overviewtable.Acquisition);
        groups = cellstr(unique(overviewtable.Combi));
        tablekw = table();
        results = table('Size',[3,2],'VariableTypes', ["string", "double"], ...
            'VariableNames', ["Acquisition", "pValue"]);
        results.Acquisition = ['A1';'A2';'A3'];
        
        for indacq = 1:size(Acquisitions, 1)
            groupdataAcq = overviewtable(overviewtable.Acquisition == Acquisitions(indacq),:);
            
            for indgroup = 1:size(groups, 1)
                group = groupdataAcq(groupdataAcq.Combi == groups(indgroup),:);
                
                if size(group,1) < size(tablekw,1)
                    warning('off','all');
                    group.Movement(end+1:size(tablekw,1)) = NaN;
                    warning('on','all');
                elseif size(group,1) > size(tablekw,1) && indgroup ~= 1
                    rowdifference = size(group,1)-size(tablekw,1);
                    tablekw{end+1:end+rowdifference,:} = missing;
                end
                
                eval(['tablekw.' groups{indgroup} ' = group.Movement;']);
            end
            p = kruskalwallis(table2array(tablekw), [],'off');
            results.pValue(indacq) = p;
        end
        
end
