% NumberOfPlots = 'All' or 'Selection'
% Side = 'L' or 'R'

function PlotUSResultsBoxplot(ResultsUS, NumberOfPlots, Side)

%% set up
% Get good tables
if matches(Side, 'R')
    ResultsUS = ResultsUS(ResultsUS.X == 1, :);
elseif matches(Side, 'L')
    ResultsUS = ResultsUS(ResultsUS.X == 1, :);
else
    disp('Side not recognized, take R as default.')
    ResultsUS = ResultsUS(ResultsUS.X == 1, :);
end

ResultsUS.xaxisplacement = ones(size(ResultsUS, 1),1);
ResultsUS(matches(ResultsUS.CaClSham, 'CaCl'), 'xaxisplacement') = 0.75.*ResultsUS(matches(ResultsUS.CaClSham, 'CaCl'), 'xaxisplacement');
ResultsUS(matches(ResultsUS.CaClSham, 'Sham'), 'xaxisplacement') = 1.25.*ResultsUS(matches(ResultsUS.CaClSham, 'Sham'), 'xaxisplacement');

% Set up titles
if matches(NumberOfPlots, 'All')
    results = {'PulsatilityIndex', 'ResistenceIndex', 'MeanVelocity',...
        'DiameterChange', 'PCTDiameterChange'};
    ylabels = {'delta velocity/mean velocity', 'delta velocity/max velocity', ...
        'mm/sec','mm','percentage'};
elseif matches(NumberOfPlots, 'Selection')
    results = {'PulsatilityIndex', 'MeanVelocity'};
    ylabels = {'delta velocity/mean velocity', 'mm/sec'};
else
    disp('NumberOfPlots not recognized, take ''All'' as default.')
    results = {'PulsatilityIndex', 'ResistenceIndex', 'MeanVelocity',...
        'DiameterChange', 'PCTDiameterChange'};
    ylabels = {'delta velocity/mean velocity', 'delta velocity/max velocity', ...
        'mm/sec','mm','percentage'};
end

% Get example frame
frames = GetFramesRecordingExample('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\US-example-data\Used - PWR1 -- 2023-03-31-15-30-26_M35 PW mode 4 right-2023-03-31-13-56-35_1.avi');

% Get pvalues
pvalues = NaN(1, size(results,2));
for ind = 1:size(results, 2)
    CaCl = table2array(ResultsUS(matches(ResultsUS.CaClSham, 'CaCl'), results{ind}));
    Sham =  table2array(ResultsUS(matches(ResultsUS.CaClSham, 'Sham'), results{ind}));

    % check normal distr. with anderson darling test and homogeneity of
    % variance with vartest2
    if adtest(CaCl) == 0 && adtest(Sham) == 0 && vartest2(CaCl, Sham) == 0
        % disp([results{ind} ' - All assumptions t-test met!']);
        [~,p] = ttest2(CaCl, Sham);
        pvalues(ind) = p;
    else
        disp([results {ind} ' - Assumptions t-test not met. Find different statistical test.'])
    end
end

%% Start plotting
f = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
if matches(NumberOfPlots, 'All')
    t = tiledlayout('flow');
elseif matches(NumberOfPlots, 'Selection')
    t = tiledlayout(3,2);
end

% Example Image
if matches(NumberOfPlots, 'All')
    nexttile
    f.Position = [10 80 1200 800];
elseif matches(NumberOfPlots, 'Selection')
    % nexttile(3,[2,2])
    nexttile(1,[2,2])
    f.Position = [10 80 800 800];
    axis off
end

imagesc(frames(:,:,2,2));
title('Example frame', 'FontSize', 15);
axis off

% Boxplots
for ind = 1:size(results,2)
    nexttile
    eval(['b = boxchart(ResultsUS.' results{ind} ', ''GroupByColor'', categorical(ResultsUS.CaClSham), ' ...
        '''LineWidth'', 2, ''MarkerStyle'', ''none'');']);
    b(1).SeriesIndex = 7;
    b(2).SeriesIndex = 1;
    hold on

    % plot individual points for mice
    eval(['hs = scatter(ResultsUS.xaxisplacement, ResultsUS.' results{ind} ', 50, ''filled'', ' ...
        '''jitter'', ''on'', ''JitterAmount'', 0.02);'])
    hs.MarkerFaceColor = [0 0 0];
    hs.MarkerFaceAlpha = 0.3;
    axes1 = b.Parent;
    hold(axes1, 'on');
    set(axes1, 'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    ylabel(ylabels{ind});
    title(results{ind}, 'FontSize', 15);
    xlabel(['p = ' num2str(pvalues(ind))]);
    xticks([]);
    xticklabels({})
end

lgd = legend(b, 'NumColumns', 2);
lgd.FontSize = 15;
lgd.Box = 'off';
lgd.Layout.Tile = 'south';

% title(t, 'Ultrasound Results Carotid Artery', 'FontSize', 20)

%% save
saveas(gcf, ['C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\USboxplot_' Side '_' NumberOfPlots '.tiff'], 'tiff');
saveas(gcf, ['C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\USboxplot_' Side '_' NumberOfPlots '.epsc'], 'epsc');

end


function frames = GetFramesRecordingExample(path_recording)

v=VideoReader(path_recording);
v.CurrentTime=1-1/30;
tmp=readFrame(v);
frames = []; %to make sure the pwv_frames is in 4D double and not uint
frames(1:size(tmp,1),:,:,1)=tmp;

if v.Duration >1 && v.Duration <2
    v.CurrentTime=1+1/30;
    tmp=readFrame(v);
    frames(1:size(tmp,1),:,:,2)=tmp;
end
if v.Duration>=2
    v.CurrentTime=2;
    tmp=readFrame(v);
    frames(1:size(tmp,1),:,:,2)=tmp;
end
if v.Duration>=3+1/30
    v.CurrentTime=3+1/30;
    tmp=readFrame(v);
    frames(1:size(tmp,1),:,:,3)=tmp;
end

end