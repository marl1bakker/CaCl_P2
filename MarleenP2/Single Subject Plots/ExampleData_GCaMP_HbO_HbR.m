%% Plot example HbO/HbR/GCaMP for fig 1
% cd 'D:\MarleenTemp\M32-A1-R2\CtxImg'
% DataFolder = 'D:\MarleenTemp\M32-A1-R2\CtxImg\';
% DataFolder = '/media/mbakker/GDrive2/P2/GCaMP/M32/A1-R2/CtxImg';

% Note: M32, acquisition 1, Recording 2 used.

function ExampleData_GCaMP_HbO_HbR(DataFolder, SaveFolder, coords, timing)
if ~matches(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist("SaveFolder", 'var')
    % SaveFolder = 'C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Article\Figures\MatlabGenerated\GCaMP\';
    SaveFolder =  '/media/mbakker/GDrive/P2/GCaMP/SingleSubjectExample/';
else
    if ~matches(SaveFolder(end), filesep)
        SaveFolder = [SaveFolder filesep];
    end
end

%% Movement
load([DataFolder 'MovMask.mat'], 'MovMask');

%% fluo data
fid = fopen([DataFolder 'hemoCorr_fluo.dat']);
dF = fread(fid, inf, '*single');
dF = reshape(dF, 512,512,[]);
fclose(fid);

if ~exist("coords", 'var')
    coords = [300 300];
end
if ~exist("timing", 'var')
    timing = [800 1100];
end

% f1 = figure;
% imagesc(dF(:,:,20))
% f1.Position = [20 30 500 400];
% f2 = figure;
% plot(squeeze(dF(coords(1),coords(2), timing(1):timing(2))))
% f2.Position = [20 450 500 400];
% f3 = figure;
% plot(MovMask)
% ylim([-0.25 1.25]);
% xline(timing(1))
% xline(timing(2))
% f3.Position = [400 100 1000 200];

%% get activation points
% fluo
%Zscore:
zF = (dF - mean(dF, 3))./std(dF,0,3);
%Threshold on Zscore:
aF = zF >= 1.95;
% Removing noise:
for ind = 1:size(aF,3)
    aF(:,:,ind) = bwmorph(bwmorph(aF(:,:,ind), 'close', inf),'open',inf);
end
%Now, we want only the beginning of activations:
aF = aF(:,:,2:end)&~aF(:,:,1:(end-1));
aF = cat(3, false(size(aF,1),size(aF,2)), aF);

clear zF ind


%% hbo data
fid = fopen([DataFolder 'HbO.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH, 512,512,[]);
fclose(fid);

% HbO spontaneous acts
zF = (dH - mean(dH, 3))./std(dH,0,3);
%Threshold on Zscore:
aF_h = zF >= 1.95;
% Removing noise:
for ind = 1:size(aF_h,3)
    aF_h(:,:,ind) = bwmorph(bwmorph(aF_h(:,:,ind), 'close', inf),'open',inf);
end
%Now, we want only the beginning of activations:
aF_h = aF_h(:,:,2:end)&~aF_h(:,:,1:(end-1));
aF_h = cat(3, false(size(aF_h,1),size(aF_h,2)), aF_h);

clear zF ind


%% hbr data 
fid = fopen([DataFolder 'HbR.dat']);
dR = fread(fid, inf, '*single');
dR = reshape(dR, 512,512,[]);
fclose(fid);


%% plot together
%one pixel
f = figure('InvertHardcopy','off', 'Color', [1 1 1]);
PlotExampleData(coords, timing, dF, aF, dH, dR, aF_h)
f.Position = [10 10 1500 800];
title('Example Data Single Pixel','FontSize', 20)

% if you want to do multiple in one figure, comment out the ylabels
% f = figure('InvertHardcopy','off', 'Color', [1 1 1]);
% t = tiledlayout(3,2);
% nexttile
% PlotExampleData([200 170], timing, dF, dH, dR)
% nexttile
% PlotExampleData([200 300], timing, dF, dH, dR)
% nexttile
% PlotExampleData([300 150], timing, dF, dH, dR)
% nexttile
% PlotExampleData([300 320], timing, dF, dH, dR)
% nexttile
% PlotExampleData([400 140], timing, dF, dH, dR)
% nexttile
% PlotExampleData([400 340], timing, dF, dH, dR)

%% Save
% saveas(gcf, [SaveFolder 'ExampleData_Pixel_' num2str(coords(1)) '-' num2str(coords(2)) '_Timing_' num2str(timing(1)) '-' num2str(timing(2)) '.tiff'], 'tiff');
% saveas(gcf, [SaveFolder 'ExampleData_Pixel_' num2str(coords(1)) '-' num2str(coords(2)) '_Timing_' num2str(timing(1)) '-' num2str(timing(2)) '.epsc'], 'epsc');
saveas(gcf, [SaveFolder 'ExampleData_Pixel_' num2str(coords(1)) '-' num2str(coords(2)) '_Timing_' num2str(timing(1)) '-' num2str(timing(2)) '.svg'], 'svg');

end


function PlotExampleData(coords, timing, dF, aF, dH, dR, aF_h)
fluopix = dF(coords(1), coords(2), timing(1):timing(2));
fluopix = reshape(fluopix, 1, []);
fluopix = fluopix - 1; %center around 0
acts = aF(coords(1), coords(2), timing(1):timing(2));
acts = reshape(acts, 1, []);
cutoff_fluo = fluopix(find(acts,1));

hbopix = dH(coords(1), coords(2), timing(1):timing(2));
hbopix = reshape(hbopix, 1, []);
acts_h = aF_h(coords(1), coords(2), timing(1):timing(2));
acts_h = reshape(acts_h, 1, []);
cutoff_hbo = hbopix(find(acts_h, 1));

hbrpix = dR(coords(1), coords(2), timing(1):timing(2));
hbrpix = reshape(hbrpix, 1, []);

hbtpix = hbopix + hbrpix;

% freq is 15 hz
nrframes = timing(2)-timing(1)+1;
x = linspace(1, nrframes/15, nrframes);
ylimfluo = [-0.05 0.1];

% Right side
yyaxis right
plot(x, hbopix, 'LineWidth', 2)
ylabel('Hemodynamics (\muM)', 'interpreter', 'Tex', 'Rotation', 270)
hold on
plot(x, hbrpix, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2)
plot(x, hbtpix, 'Color', 'black', 'LineStyle', '-', 'LineWidth', 2)
% ylim([-4.29 8.56]); % these are weird values to have the 0 of hemo at the same level as 1 of the fluo
ylim(ylimfluo*85)
yticks([-4 -2 0 2 4 6 8])
ax = gca;
ax.YColor = 'red';
ax.XColor = 'k';
% threshold
line([x(1) x(end)], [cutoff_hbo cutoff_hbo], 'LineWidth' , 1, 'LineStyle', '--', 'Color', [0.6350 0.0780 0.1840])

%Left side
yyaxis left
plot(x, fluopix, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
ylabel('GCaMP (\Delta F/F)')
ylim(ylimfluo)
ax = gca;
% yticks([0.97 1.00 1.03 1.06 1.09])
yticks([-0.05 -0.025 0 0.025 0.05 0.075 0.1])
ax.YColor = [0.4660 0.6740 0.1880];
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)
% threshold
hold on
line([x(1) x(end)], [cutoff_fluo cutoff_fluo], 'LineWidth' , 1, 'LineStyle', '--', 'Color', [0.4660 0.6740 0.1880]);

xlabel('Time (seconds)')
legend({'GCaMP', 'GCaMP threshold', 'HbO', 'HbR', 'HbT', 'HbO threshold'})
xlim([0.5 20.5])
subtitle(['Timing is ' num2str(timing(1)) ' - ' num2str(timing(2)) ' frames, pixel is Y: ' num2str(coords(1)) ', X: ' num2str(coords(2))]);
end