function DVARS(DataFolder)

SaveFolder = '/media/mbakker/GDrive/P2/GCaMP/Movement/';

%% set up
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end


%% get data
datatypes = {'hemoCorr_fluo', 'HbO', 'HbR'};
dvars = [];

for inddata = 1:length(datatypes)
    fid = fopen([DataFolder datatypes{inddata} '.dat']);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, 512*512, []);
    fclose(fid);

    % each frame minus frame before
    dvars_data = NaN(size(dat));

    for ind = 2:size(dat,2)
        dvars_data(:,ind) = dat(:,ind)-dat(:,ind-1);
    end

    % square pixels
    dvars_data = dvars_data.^2;

    % average per frame
    dvars_data = mean(dvars_data, 1, 'omitnan');

    % square root
    dvars_data = sqrt(dvars_data);

    dvars.(datatypes{inddata}) = dvars_data;

    disp(['... ' datatypes{inddata}, ' done ...'])
end

%% plot all
load([DataFolder 'MovMask.mat'], 'MovMask')
load([DataFolder 'OutlierMask.mat'], 'OutlierFrames')

f = figure;
plot(MovMask.*OutlierFrames, 'Color', 'black')
hold on

plot(dvars.hemoCorr_fluo.*10, 'Color', 'green', 'LineStyle','-') % very low, do times 10
plot(dvars.HbO, 'Color', 'red', 'LineStyle','-')
plot(dvars.HbR, 'Color', 'blue', 'LineStyle','-')

legend({'Movmask+OutlierFrames', 'GCaMP DVARS * 10', 'HbO DVARS', 'HbR DVARS'})

seps = strfind(DataFolder, filesep);
Mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);
Acq = DataFolder(seps(end-2)+1:seps(end-1)-1);
title([Mouse ' ' Acq])

f.Position = [200 400 1500 400];

saveas(f, [SaveFolder 'DVARS_' Mouse '_' Acq '.png'], 'png')
saveas(f, [SaveFolder 'DVARS_' Mouse '_' Acq '.fig'], 'fig')
close(f)
end

