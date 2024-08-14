function ExampleData_NoFilt(DataFolder, SaveDir)

% set up
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = strcat(DataFolder, filesep);
end

seps = strfind(DataFolder, filesep);
Mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);

if ~exist('SaveDir', 'var')
    SaveDir =  '/media/mbakker/GDrive/P2/GCaMP';
end

if ~exist([SaveDir '/SingleSubjectExample/'], 'dir')
    mkdir([SaveDir '/SingleSubjectExample/'])
end



%% make none filtered hemodynamics data and get average curves for those
if ~exist([DataFolder 'HbO_nofilt.dat'], 'file')
    HemoCompute_nofilt(DataFolder, DataFolder, 'gcamp', {'red','green'},1);
    GetAverageCurves(DataFolder, 'fluo_nofilt', 'NVC_NoneFilt', 'hemoCorr_fluo', 'HbO_nofilt', 'HbR_nofilt')
end

%% load both NVC files
load([DataFolder 'NVC_ROI.mat'], 'fluocurves', 'hbocurves', 'hbrcurves');
fluo.filt = fluocurves;
hbo.filt = hbocurves;
hbr.filt = hbrcurves;
% hbt.filt = hbtcurves;

load([DataFolder 'NVC_NoneFilt.mat'], 'fluocurves', 'hbocurves', 'hbrcurves');
fluo.nofilt = fluocurves;
hbo.nofilt = hbocurves;
hbr.nofilt = hbrcurves;
% hbt.nofilt = hbtcurves;

clear fluocurves hbocurves hbrcurves hbtcurves

%% Plot whole brain
f = figure('InvertHardcopy','off','Color',[1 1 1]);
x = linspace(-5, 10, 226);
opacity_filt = 1;
opacity_nofilt = 0.3;

% hbo
yyaxis right
plot(x,hbo.filt.WholeBrain, 'Color', [1 0 0 opacity_filt], 'LineStyle', '-', 'LineWidth', 2)
hold on
plot(x,hbo.nofilt.WholeBrain, 'Color', [1 0 0 opacity_nofilt], 'LineStyle', '-', 'LineWidth', 2)

% hbr
plot(x,hbr.filt.WholeBrain, 'Color', [0 0 1 opacity_filt], 'LineStyle', '-', 'LineWidth', 2)
plot(x,hbr.nofilt.WholeBrain, 'Color', [0 0 1 opacity_nofilt], 'LineStyle', '-', 'LineWidth', 2)

h_label = ylabel('Hemodynamics (DmM)','Rotation', 270); % to work with illustrator
% h_label = ylabel('Hemodynamics (\Delta \muM)', 'interpreter', 'Tex','Rotation', 270);

ax = gca;
ax.YColor = 'red';
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
xlabel('Time (sec.)')
ylim([-2 2])

% fluo
yyaxis left
plot(x,fluo.filt.WholeBrain-1, 'Color', [0.4660 0.6740 0.1880 opacity_filt], 'LineStyle','-','LineWidth',2)
plot(x,fluo.nofilt.WholeBrain-1, 'Color', [0.4660 0.6740 0.1880 opacity_nofilt], 'LineStyle','-','LineWidth',2)

f_label = ylabel('GCaMP (DF/F)'); % to work with illustrator
% f_label = ylabel('GCaMP (\Delta F/F)', 'interpreter', 'Tex'); 

ylim([-0.05 0.05]);
ax = gca;
ax.YColor = [0.4660 0.6740 0.1880];
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)

legend({'Filtered GCaMP', 'Unfiltered GCaMP', 'Filtered HbO', 'Unfiltered HbO',...
    'Filtered HbR', 'Unfiltered HbR'})

title('Whole Brain - Single Mouse')
f.Position = [20 20 900 800];
pause(0.05)

saveas(gcf, [SaveDir '/SingleSubjectExample/Filt_vs_Unfilt_WholeBrain_' Mouse '.svg'], 'svg')
close(f)
end