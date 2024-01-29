%% Analysis in the app
% First, do us_analysis via the GUI. Save the results in a CSV file. Import
% the CSV file into matlab, and call it "ResultsUS". This is important to
% make the table later. 

%% Make a nice table
load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS.mat');

load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\Mice.mat')

ResultsUS.Mouse = cellstr(repmat('empty', size(ResultsUS,1),1));
ResultsUS.Side = cellstr(repmat('empty', size(ResultsUS,1), 1));
ResultsUS.CaClSham = cellstr(repmat('empty', size(ResultsUS,1),1));
ResultsUS.MaleFemale = cellstr(repmat('empty', size(ResultsUS,1),1));

for ind = 1:size(ResultsUS, 1)
    fullname = ResultsUS.Row{ind};
    ResultsUS.Mouse(ind,:) = cellstr(fullname(1:strfind(fullname, '-')-1));
    ResultsUS.Side(ind,:) = cellstr(fullname(strfind(fullname, '-')+1:end));

    indmouse = find(ismember(Mice.CodeOfMouse, ResultsUS.Mouse(ind))); %find mouse
    ResultsUS.CaClSham(ind,:) = cellstr(Mice.CaClSham(indmouse));
    ResultsUS.MaleFemale(ind,:) = cellstr(Mice.MaleFemale(indmouse));
end

clear ind indmouse fullname Mice

ResultsUS.X = single(ismember(ResultsUS.Side, 'R')); %L = 0, R = 1

save('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS.mat');

%% Analysis "by hand"
% For acquisitions that cannot be neatly analysed via the GUI, go to the 
% file US_analysis_by_hand. Pick out the right recordings and add them
% manually to the table. 
load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS+byhand.mat', 'ResultsUS');

%% Plot the results. 
PlotUSResults(ResultsUS)

saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\UltrasoundPlots.tiff', 'tiff');
saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\UltrasoundPlots.epsc', 'epsc');






%% Statistics
pvaluesRight = [];
pvaluesLeft = [];
%% PW Right
Right = ResultsUS(ResultsUS.X == 1, :);
results = {'PulsatilityIndex', 'ResistenceIndex', 'MeanVelocity'};

for ind = 1:size(results, 2)
    CaClR = table2array(Right(matches(Right.CaClSham, 'CaCl'), results{ind}));
    ShamR =  table2array(Right(matches(Right.CaClSham, 'Sham'), results{ind}));

    % check normal distr. with anderson darling test and homogeneity of
    % variance with vartest2
    if adtest(CaClR) == 0 && adtest(ShamR) == 0 && vartest2(CaClR, ShamR) == 0
        % disp([results{ind} ' - All assumptions t-test met!']);
        [h,p] = ttest2(CaClR, ShamR);
        % disp(['p = ' num2str(p) 'so h = ' num2str(h)])
        pvaluesRight(end+1) = p;
    else
        disp([results {ind} ' - Assumptions t-test not met. Find different statistical test.'])
    end
end

%% MMode Right
results = {'DiameterChange', 'PCTDiameterChange'};

for ind = 1:size(results, 2)
    CaClR = table2array(Right(matches(Right.CaClSham, 'CaCl'), results{ind}));
    ShamR =  table2array(Right(matches(Right.CaClSham, 'Sham'), results{ind}));

    % check normal distr. with anderson darling test and homogeneity of
    % variance with vartest2
    if adtest(CaClR) == 0 && adtest(ShamR) == 0 && vartest2(CaClR, ShamR) == 0
        % disp([results{ind} ' - All assumptions t-test met!']);
        [h,p] = ttest2(CaClR, ShamR);
        % disp(['p = ' num2str(p) 'so h = ' num2str(h)])
        pvaluesRight(end+1) = p;
    else
        disp([results {ind} ' - Assumptions t-test not met. Find different statistical test.'])
    end
end


%% plot only right side, boxplot
Right = ResultsUS(ResultsUS.X == 1, :);
Right.xaxisplacement = ones(size(Right, 1),1);
Right(matches(Right.CaClSham, 'CaCl'), 'xaxisplacement') = 0.75.*Right(matches(Right.CaClSham, 'CaCl'), 'xaxisplacement');
Right(matches(Right.CaClSham, 'Sham'), 'xaxisplacement') = 1.25.*Right(matches(Right.CaClSham, 'Sham'), 'xaxisplacement');

results = {'PulsatilityIndex', 'ResistenceIndex', 'MeanVelocity',...
    'DiameterChange', 'PCTDiameterChange'};
ylabels = {'delta velocity/mean velocity', 'delta velocity/max velocity', ...
    'mm/sec','mm','percentage'};


f = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
t = tiledlayout(2,3);
f.Position = [10 80 1200 600];

for ind = 1:size(results,2)
    nexttile
    % eval(['b = boxchart(categorical(Right.CaClSham), Right.' results{ind} ', ' ...
    %     '''LineWidth'', 2, ''MarkerStyle'', ''none'');'])
    eval(['b = boxchart(Right.' results{ind} ', ''GroupByColor'', categorical(Right.CaClSham), ' ...
        '''LineWidth'', 2, ''MarkerStyle'', ''none'');']);
    b(1).SeriesIndex = 7;
    b(2).SeriesIndex = 1;
    hold on
    eval(['hs = scatter(Right.xaxisplacement, Right.' results{ind} ', 50, ''filled'', ' ...
        '''jitter'', ''on'', ''JitterAmount'', 0.02);'])
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.3;
    % xlim([0.7 1.3]);
    axes1 = b.Parent;
    hold(axes1, 'on');
    set(axes1, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    ylabel(ylabels{ind});
    title(results{ind}, 'FontSize', 15);
    xlabel(['p = ' num2str(pvaluesRight(ind))]);
    xticks([]);
    xticklabels({})
end

lgd = legend([b]);
lgd.Layout.Tile = 6;
lgd.FontSize = 15;
lgd.Box = 'off';
% lgd.Layout.Tile = 'east'

title(t, 'Ultrasound Results Right Carotid Artery', 'FontSize', 20)

saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\USboxplotRight.tiff', 'tiff');
saveas(gcf, 'C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\USboxplotRight.epsc', 'epsc');

%% PW Left
Left = ResultsUS(ResultsUS.X == 0, :);
results = {'PulsatilityIndex', 'ResistenceIndex', 'MeanVelocity'};

for ind = 1:size(results, 2)
    CaClR = table2array(Left(matches(Left.CaClSham, 'CaCl'), results{ind}));
    ShamR =  table2array(Left(matches(Left.CaClSham, 'Sham'), results{ind}));

    % check normal distr. with anderson darling test and homogeneity of
    % variance with vartest2
    if vartest2(CaClR, ShamR) == 0
        disp([results{ind} ' - Left - All assumptions t-test met!']);
        [h,p] = ttest2(CaClR, ShamR);
        disp(['p = ' num2str(p) 'so h = ' num2str(h)])
        pvaluesLeft(end+1) = p;
    else
        disp([results {ind} ' - Assumptions t-test not met. Find different statistical test.'])
    end
end

%% MMode Left
results = {'DiameterChange', 'PCTDiameterChange'};

for ind = 1:size(results, 2)
    CaClR = table2array(Left(matches(Left.CaClSham, 'CaCl'), results{ind}));
    ShamR =  table2array(Left(matches(Left.CaClSham, 'Sham'), results{ind}));

    % check normal distr. with anderson darling test and homogeneity of
    % variance with vartest2
    if vartest2(CaClR, ShamR) == 0
        disp([results{ind} ' - Left - All assumptions t-test met!']);
        [h,p] = ttest2(CaClR, ShamR);
        disp(['p = ' num2str(p) 'so h = ' num2str(h)])
        pvaluesLeft(end+1) = p;
    else
        disp([results {ind} ' - Assumptions t-test not met. Find different statistical test.'])
    end
end


