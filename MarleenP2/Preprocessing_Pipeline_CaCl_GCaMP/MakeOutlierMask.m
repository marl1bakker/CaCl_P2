% Takes outlier FRAMES and makes a mask for them. It calculates this
% without removing the movement mask first, so there will very likely be an
% overlap. It takes M26, A2 seperately since this acquisition has a weird
% artefact.

function MakeOutlierMask(DataFolder, datatypes, overwrite)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datatypes', 'var')
    datatypes = {'hemoCorr_fluo', 'HbO', 'HbR'};
end

if ~exist('overwrite', 'var')
    overwrite = 0;
end

% if ~exist('manualinput', 'var')
%     manualinput = 0;
% end

%% get masks
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;
load([DataFolder 'MovMask.mat'], 'MovMask');
Mouse = DataFolder(seps(end-3)+1:seps(end-2)-1);
Acq = DataFolder(seps(end-2)+1:seps(end-2)+2);

if exist([DataFolder 'OutlierMask.mat']) 
    infovar = who('-file', [DataFolder 'OutlierMask.mat']);
else
    infovar = 'none';
end

%% per frame
if contains(infovar, 'OutlierMask') % to correct for old, can delete later
    load([DataFolder 'OutlierMask.mat'], 'OutlierMask')
    OutlierFrames = OutlierMask;
    clear OutlierMask
    save([DataFolder 'OutlierMask.mat'], 'OutlierFrames');
elseif  sum(contains(infovar, 'OutlierFrames')) && overwrite == 0
    load([DataFolder 'OutlierMask.mat'], 'OutlierFrames')
else
    progress = 0;
    w = waitbar(progress, 'Outlier Masks Frames...');
    OutlierFrames = zeros(size(MovMask));

    for ind = 1:size(datatypes,2)
        datatype = datatypes{ind};

        % load data
        fid = fopen([DataFolder datatype '.dat']);
        dat = fread(fid, inf, '*single');
        dat = reshape(dat, 512,512,[]);
        fclose(fid);

        %take only brain
        dat = dat.*mask;
        dat(dat == 0) = NaN;
        dat = reshape(dat,[], size(dat,3));

        % Will take outliers based on mean, take M26 seperately (light artefact)
        % Checked based on median also, but it takes many frames which seemed
        % to be good.
        meandat = mean(dat,1,'omitnan'); % get average per frame
        indoutliers = isoutlier(meandat, "mean", 2);

        if matches(Mouse, 'M26') && matches(Acq, 'A2')
            indoutliers(7400:8700) = 1;
        end

        OutlierFrames = OutlierFrames + indoutliers;

        progress = progress + 0.3;
        waitbar(progress, w, 'Outlier Mask Frames...');
    end
    close(w)

    OutlierFrames(OutlierFrames>0) = 1;
    OutlierFrames=~OutlierFrames; %flip, 0 is outlier (should be removed, so can do dat.*mask)

    save([DataFolder 'OutlierMask.mat'], 'OutlierFrames');

end

%% per pixel
% have to do this after doing outlierframes for all datatypes, because
% outlierframes are usually because of movement or weird light artefacts
% that influence all channels, but pixels or window artefacts can be
% specific to certain lights. 

%% patch - delete later!
% this is just to keep the old outliermask should the new one give problems
% load([DataFolder 'OutlierMask.mat'], 'OutlierPixels_old');
% if ~exist('OutlierPixels_old', 'var')
%     load([DataFolder 'OutlierMask.mat'], 'OutlierPixels')
%     OutlierPixels_old = OutlierPixels;
%     save([DataFolder 'OutlierMask.mat'], 'OutlierPixels_old', '-append');
% end

%%
if sum(contains(infovar, 'OutlierPixels')) && overwrite == 0
    disp('Already done')
    return
else
    progress = 0;
    w = waitbar(progress, 'Outlier Masks Pixels...');

    for ind = 1:size(datatypes,2)
        datatype = datatypes{ind};

        % load data
        fid = fopen([DataFolder datatype '.dat']);
        dat = fread(fid, inf, '*single');
        dat = reshape(dat, 512,512,[]);
        fclose(fid);

        %% per pixel
        %take only brain
        dat = dat.*mask;
        %get rid of outlier frames
        dat = reshape(dat, 512*512, []);
        dat = dat.*MovMask;
        dat = dat.*OutlierFrames;
        dat(dat == 0) = NaN;
        dat = reshape(dat, 512, 512, []);

        %% new 11/6/24
        outlierpix = sum(isoutlier(dat, 'median', 1),3) + ...
            sum(isoutlier(dat, 'median', 2),3);
        threshold = 1000;
        OutlierPixels.(datatype) = outlierpix<threshold; %keep good ones
        OutlierPixels.raw.(datatype) = outlierpix;
        OutlierPixels.thresholds.(datatype) = threshold;
        f = figure;
        imagesc(outlierpix)
        f2 = figure;
        imagesc(outlierpix<threshold)
        close(f, f2)

        %% old
        % % get outlier pixels
        % % outlierpix = sum(isoutlier(dat, 'median',3),3);
        % outlierpix = sum(isoutlier(dat, 'median', 3),3);
        % threshold = 300;
        % 
        % % check
        % f1 = figure;
        % % imagesc(sum(outlierpix, 3), [0 300])
        % imagesc(outlierpix, [0 300]);
        % title('outliersum')
        % f1.Position = [20 20 500 400];
        % f2 = figure;
        % imagesc(dat(:,:,find(MovMask, 1, 'first')))
        % title('dat single frame')
        % f2.Position = [520 20 500 400];
        % f3 = figure;
        % imagesc(mean(dat,3, 'omitnan'))
        % title('average dat')
        % f3.Position = [1020 20 500 400];
        % 
        % outlierpix(outlierpix<threshold) = 1;
        % outlierpix(outlierpix>=threshold) = 0;
        % f4 = figure;
        % imagesc(outlierpix)
        % title('Outlier mask')
        % f4.Position = [1020 450 500 400];
        % 
        % % answer = questdlg(['Outliermask makes sense? ' datatype]);
        % % while matches(answer, 'No')
        % %     % make new threshold, check if it's better
        % %     outlierpix = sum(isoutlier(dat, 'median', 3),3);
        % %     threshold = inputdlg('New Threshold: ');
        % %     threshold = str2double(threshold);
        % % 
        % %     close(f4)
        % %     outlierpix(outlierpix<threshold) = 1;
        % %     outlierpix(outlierpix>=threshold) = 0;
        % %     f4 = figure;
        % %     imagesc(outlierpix)
        % %     title('Outlier mask')
        % %     f4.Position = [1020 450 500 400];
        % % 
        % %     answer = questdlg(['Outliermask makes sense? ' datatype]);
        % %     % return
        % % end
        % % if matches(answer, 'Cancel')
        % %     return
        % % end
        % 
        % close(f1, f2, f3, f4)
        % eval(['OutlierPixels_old.' datatype ' = outlierpix;'])
        % eval(['OutlierPixels_old.thresholds.' datatype ' = threshold;'])
        % 
        % 
        % % answer = questdlg('Manual adjustments needed?');
        % % 
        % % while matches(answer, 'Yes')
        % %     imagesc(mean(dat,3, 'omitnan'))
        % %     addroi = drawpolygon;
        % % 
        % %     bubbles = bubbles + addroi;
        % % end
        % % 
        % % eval(['ManualOutliers.' datatype ' = ;'])

        progress = progress + 0.15;
        waitbar(progress, w, 'Outlier Mask Pixels...');


    end
end


%% save
close(w)
% save([DataFolder 'OutlierMask.mat'], 'OutlierPixels_old', '-append');
save([DataFolder 'OutlierMask.mat'], 'OutlierPixels', '-append');

end