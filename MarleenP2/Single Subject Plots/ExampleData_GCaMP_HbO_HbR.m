%% Plot example HbO/HbR/GCaMP for fig 1
% cd 'D:\MarleenTemp\M32-A1-R2\CtxImg'
% DataFolder = 'D:\MarleenTemp\M32-A1-R2\CtxImg\';

% Note: M32, acquisition 1, Recording 2 used.

function ExampleData_GCaMP_HbO_HbR(DataFolder, SaveFolder, coords, timing)
if ~matches(DataFolder(end), filesep)
    DataFolder = [DataFolder filesep];
end

if ~exist("SaveFolder", 'var')
    SaveFolder = 'C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Article\Figures\MatlabGenerated\GCaMP\';
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

% for ind = 1:500
%     imagesc(dF(:,:,ind), [0.95 1.05])
%     title(num2str(ind))
%     pause(0.05)
% end

if ~exist("coords", 'var')
    coords = [300 300];
end
if ~exist("timing", 'var')
    timing = [800 1100];
end

f1 = figure;
imagesc(dF(:,:,20))
f1.Position = [20 30 500 400];
f2 = figure;
plot(squeeze(dF(coords(1),coords(2), timing(1):timing(2))))
f2.Position = [20 450 500 400];
f3 = figure;
plot(MovMask)
ylim([-0.25 1.25]);
xline(timing(1))
xline(timing(2))
f3.Position = [400 100 1000 200];

% fluopix = dF(coords(1), coords(2), timing(1):timing(2));
% fluopix = reshape(fluopix, 1, []);

%% hbo data
fid = fopen([DataFolder 'HbO.dat']);
dH = fread(fid, inf, '*single');
dH = reshape(dH, 512,512,[]);
fclose(fid);

% for ind = 400:600
%     imagesc(dH(:,:,ind), [-40 40])
%     title(num2str(ind))
%     pause(0.05)
% end

%% hbr data 
fid = fopen([DataFolder 'HbR.dat']);
dR = fread(fid, inf, '*single');
dR = reshape(dR, 512,512,[]);
fclose(fid);

% for ind = 400:600
%     imagesc(dR(:,:,ind), [-10 10])
%     title(num2str(ind))
%     pause(0.05)
% end


%% plot together
%one pixel
f = figure('InvertHardcopy','off', 'Color', [1 1 1]);
PlotExampleData(coords, timing, dF, dH, dR)
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
saveas(gcf, [SaveFolder 'ExampleData_Pixel_' num2str(coords(1)) '-' num2str(coords(2)) '_Timing_' num2str(timing(1)) '-' num2str(timing(2)) '.tiff'], 'tiff');
saveas(gcf, [SaveFolder 'ExampleData_Pixel_' num2str(coords(1)) '-' num2str(coords(2)) '_Timing_' num2str(timing(1)) '-' num2str(timing(2)) '.epsc'], 'epsc');

end


function PlotExampleData(coords, timing, dF, dH, dR)
fluopix = dF(coords(1), coords(2), timing(1):timing(2));
fluopix = reshape(fluopix, 1, []);
hbopix = dH(coords(1), coords(2), timing(1):timing(2));
hbopix = reshape(hbopix, 1, []);
hbrpix = dR(coords(1), coords(2), timing(1):timing(2));
hbrpix = reshape(hbrpix, 1, []);

% freq is 15 hz
nrframes = timing(2)-timing(1)+1;
x = linspace(1, nrframes/15, nrframes);

%Left side
yyaxis left
plot(x, fluopix, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
ylabel('Change in GCaMP fluorecence (\Delta F/F)')
ylim([0.95 1.1]);
yticks([0.97 1.00 1.03 1.06 1.09])
ax = gca;
ax.YColor = 'k';
ax.XColor = 'k';
set(ax, 'FontSize', 15, 'LineWidth', 2)

% Right side
yyaxis right
plot(x, hbopix, 'LineWidth', 2)
ylabel('Change in hemodynamics (\muM)', 'interpreter', 'Tex', 'Rotation', 270)
hold on
plot(x, hbrpix, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2)
ylim([-4.29 8.56]); % these are weird values to have the 0 of hemo at the same level as 1 of the fluo
yticks([-2 0 2 4 6 8])
ax.YColor = 'k';

xlabel('Time (seconds)')
legend({'GCaMP', 'HbO', 'HbR'})
xlim([0.5 20.5])
subtitle(['Timing is ' num2str(timing(1)) ' - ' num2str(timing(2)) ' frames, pixel is Y: ' num2str(coords(1)) ', X: ' num2str(coords(2))]);
end