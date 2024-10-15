function CorrMatricesStatsA1A2A3(GSR)

if ~exist('GSR', 'var') || GSR == 0
    % GSR = 0;
    addgsr = '';
elseif GSR == 1
    addgsr = '_GSR';
end

Groups = {'CaCl', 'Sham'};
SaveFolder = '/media/mbakker/GDrive/P2/GCaMP/CombinedCorrMatrices/';
datatypes = {'GCaMP', 'HbO', 'HbR'};
Acquisitions = {'A1', 'A2', 'A3'};
allrois = {'VR', 'SR', 'MR', 'RR', 'VL', 'SL', 'ML', 'RL'};
[allrois, sortindex] = sort(allrois); % alfabetisch

f = figure('InvertHardcopy','off','Color',[1 1 1]);
t = tiledlayout(1, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

% go per datatype
for inddata = 1:size(datatypes, 2)
    datatype = datatypes{inddata};

    %% get data
    % go per group
    for indgroup = 1:size(Groups, 2)
        group = Groups{indgroup};

        % get the corrmatrices of all acq of one group
        for indAcq = 1:size(Acquisitions, 2)
            Acq = Acquisitions{indAcq};

            % Get data
            load([SaveFolder datatype filesep group '_' Acq '_corrmatrices' addgsr '.mat'], 'allcorrmatrices')

            % z transform fisher
            allcorrmatrices = atanh(allcorrmatrices);

            allcorrmatrices(allcorrmatrices == inf) = NaN; %mouse x roi
            matrices.(group).(Acq) = allcorrmatrices;
        end

    end


    %% stats
    % Compare for one datatype (e.g. GCaMP), two group (cacl sham), if
    % there is a sign. effect of time. This has to be done per seed
    % pair.

    pvalue_time = NaN(size(allrois,2)^2, 1);
    % pvalue_group_time = NaN(size(allrois,2)^2, 1);

    for indseedpair = 1:size(allrois,2)^2
        % if you only have nan because it's the seed with itself, skip
        if ~any(matrices.CaCl.A1(:,indseedpair))
            continue
        end

        % build good table
        cacltable = table('Size', [size(matrices.CaCl.A1,1) 4], ...
            'VariableNames', {'Group', 'A1', 'A2', 'A3'}, ...
            'VariableTypes', {'categorical', 'double', 'double','double'});
        cacltable.Group = repmat('CaCl', size(cacltable,1),1);
        cacltable.A1 = matrices.CaCl.A1(:,indseedpair);
        cacltable.A2 = matrices.CaCl.A2(:,indseedpair);
        cacltable.A3 = matrices.CaCl.A3(:,indseedpair);

        shamtable = table('Size', [size(matrices.Sham.A1,1) 4], ...
            'VariableNames', {'Group', 'A1', 'A2', 'A3'}, ...
            'VariableTypes', {'categorical', 'double', 'double','double'});
        shamtable.Group = repmat('Sham', size(shamtable,1),1);
        shamtable.A1 = matrices.Sham.A1(:,indseedpair);
        shamtable.A2 = matrices.Sham.A2(:,indseedpair);
        shamtable.A3 = matrices.Sham.A3(:,indseedpair);

        seedpairtable = [cacltable; shamtable]; % table with all connectivities of a single seedpair
        clear shamtable cacltable


        %% check assumptions ranova
        % check normal distribution (kolmogorov smirnov)
        [h1, ~] = kstest(seedpairtable.A1);
        [h2, ~] = kstest(seedpairtable.A2);
        [h3, ~] = kstest(seedpairtable.A2);

        if h1+h2+h3 > 0
            % disp('not normal distr. Do Friedman')
            dofriedman = 1;
        else
            dofriedman = 0;
        end

        % check homogeneity of variance (bartlett)
        [p1] = vartestn(seedpairtable.A1, 'Display', 'off');
        [p2] = vartestn(seedpairtable.A2, 'Display', 'off');
        [p3] = vartestn(seedpairtable.A3, 'Display', 'off');
        if p1 < 0.05 || p2 < 0.05 || p3 < 0.05
            dofriedman = 1;
        end
        clear h1 h2 h3 p1 p2 p3



        %% start friedman
        if dofriedman

            if anynan(table2array(seedpairtable(:,2:end)))
                nanind = find(isnan(seedpairtable.A1)+isnan(seedpairtable.A2)+isnan(seedpairtable.A3)); %find the nan
                seedpairtable(nanind,:) = []; %get rid of mice with nan
            end

            pvalue_time(indseedpair) = friedman(table2array(seedpairtable(:,2:end)), 1, 'off'); % 1 repetition

        else

            error('not friedman - think about what to do')
            %% Do repeated measures anova
            Time = [1 2 3]';
            rm = fitrm(seedpairtable, 'A1-A3 ~ Group', 'WithinDesign',Time);
            ranovatbl = ranova(rm);

            pvalue_time(indseedpair) = ranovatbl.pValue(1);
            % pvalue_group_time(indseedpair) = ranovatbl.pValue(2);
        end
    end

    %% Do FDR
    pvalue_time = reshape(pvalue_time, 8, 8);
    pvalue_time = pvalue_time(sortindex, sortindex); % this sorts it in the right manner (ml mr rl rr etc)
    pvalue_time = tril(pvalue_time, -1);
    pvalue_time = reshape(pvalue_time, 64, 1);

    pvalues28 = [];
    for ind = 1:size(pvalue_time,1)
        if pvalue_time(ind) ~= 0
            pvalues28 = [pvalues28; pvalue_time(ind)];
        end
    end
    clear ind

    qvalues = mafdr(reshape(pvalues28, [], 1),'BHFDR', 'true');

    q = tril(ones(8), -1);
    q = reshape(q, 64, 1);
    ind2 = 1;
    for ind = 1:64
        if q(ind) ~= 0
            q(ind) = qvalues(ind2);
            ind2 = ind2 + 1;
        end
    end
    q = reshape(q, 8, 8);
    pvalue_time = reshape(pvalue_time, 8,8);
    pqcombi = pvalue_time' + q;
    q(q==0) = NaN;
    pvalue_time(pvalue_time==0) = NaN;

    if sum(any(q<0.05))
        disp('HHHHHHHAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA HOLY SHIT')
        % warning('YOU HAVE SIGNIFICANCE')
    end


    %% plot
    nexttile
    imagesc(pqcombi, [0 0.05])
    colormap('gray')
    title(datatype)
    axis('image')
    axis('equal')
    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca, 'YTickLabel', ay);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)

    %% save p/q values
    pqcombi(pqcombi == 0) = 1;
    pqcombi = [['ROIs'; allrois'], [allrois; num2cell(pqcombi)]];
    writecell(pqcombi, [SaveFolder 'pqvalues_A1A2A3' addgsr '.xlsx'], 'Sheet', [datatype]);

    %% make ready for next
    clear Acq allcorrmatrices datatype dofriedman group ind ind2 indAcq inddata indgroup indseedpair matrices nanind pqcombi pvalue_group_time pvalue_time pvalues28 q qvalues seedpairtable

end

close(f)

end
