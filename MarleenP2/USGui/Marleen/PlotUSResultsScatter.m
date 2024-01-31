% ResultsUS is table
% yVariable is 'DiameterChange', 'PulsatilityIndex', etc.
function PlotUSResultsScatter(ResultsUS, SaveDir)
if ~exist("SaveDir", "var")
    SaveDir = 'C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Article\Figures\MatlabGenerated\US';
end

f = figure('InvertHardcopy', 'off', 'Color', [1 1 1]);
t = tiledlayout(2,3);
f.Position = [10 10 1200 600];

% PW mode (3)
nexttile
MakeSingleUSPlot(ResultsUS, 'PulsatilityIndex');
ylabel('delta velocity/mean velocity')
nexttile
MakeSingleUSPlot(ResultsUS, 'ResistenceIndex');
ylabel('delta velocity/max velocity')
nexttile
MakeSingleUSPlot(ResultsUS, 'MeanVelocity');
ylabel('mm/sec')

% M Mode (2)
nexttile
MakeSingleUSPlot(ResultsUS, 'DiameterChange');
ylabel('mm')
nexttile
MakeSingleUSPlot(ResultsUS, 'PCTDiameterChange');
ylabel('percentage')

title(t, 'Ultrasound Results', 'FontSize', 20)

%%
leg = legend('CaCl', 'Sham');
leg.Layout.Tile = 6;

%% save
saveas(gcf, [SaveDir '\USScatterplot_BothSides.tiff'], 'tiff');
saveas(gcf, [SaveDir '\USScatterplot_BothSides.epsc'], 'epsc');

end

function MakeSingleUSPlot(ResultsUS, yVariable)
hold on
Mice = unique(ResultsUS.Mouse);

for m = 1:size(Mice,1)
    mouse = Mice{m};
    mousetable = ResultsUS(find(ismember(ResultsUS.Mouse, mouse)),:);
    jitteramount = rand(1,1)*0.1-0.05;

    eval(['p = plot(mousetable.X + jitteramount, mousetable.' yVariable ');'])
    p.LineWidth = 1;
    p.LineStyle = "--";
    if matches(mousetable.CaClSham{1}, 'CaCl')
        % p.Color = 'red';
        p.Color = [0.6350 0.0780 0.1840];
        p.Marker = '.';
        p.MarkerSize = 20;
    else
        p.Color = [0.8500 0.3250 0.0980]; %orange, for NaCl, to fit with other plots
        p.Marker = "o";
        p.MarkerSize = 4;
        % p.Color = [0.3010 0.7450 0.9330]; %light blue
        % p.Color = 'blue';
    end

    xlim([-0.3 1.3]);
    xticks([0 1]);
    xticklabels({'Left', 'Right'});

end

axes1 = p.Parent;
hold(axes1, "on")
set(axes1, 'FontSize', 12, 'LineWidth', 2);
title(yVariable, 'FontSize', 15)
hold off

end
