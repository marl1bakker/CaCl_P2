% Note: first part is same as GetAverageCurves. Would be smart to integrate
% these two.

function Brainmap_specs(DataFolder)
    % , NVCname, datanamefluocurves, datanamehbo, datanamehbr, datanamefluoactivations)


%% set up
if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('datanamefluoactivations', 'var')
    datanamefluoactivations = 'hemoCorr_fluo';
end

if ~exist('datanamefluocurves', 'var')
    datanamefluocurves = 'hemoCorr_fluo';
end

% if ~exist('NVCname', 'var') && matches(datanamefluocurves, 'hemoCorr_fluo')
%     NVCname = 'NVC_ROI';
% elseif ~exist('NVCname', 'var')
%     NVCname = ['NVC_ROI_' datanamefluocurves];
% end

if ~exist('datanamehbo', 'var')
    datanamehbo = 'HbO';
    datanamehbr = 'HbR';
end

if ~exist([DataFolder datanamefluoactivations '.dat'], 'file')
    disp([datanamefluoactivations ' could not be found, function exited'])
    return
end

if ~exist([DataFolder datanamehbo '.dat'], 'file')
    disp([datanamehbo ' could not be found, function exited'])
    return
end


%% timewindow
imfreq = 15;
timebefore = 5*imfreq; %5 sec before
timeafter = 10*imfreq; % 10 sec after
windowtime = timebefore+timeafter+1;

%% get points of activation
% load data
fid = fopen([DataFolder datanamefluoactivations '.dat']);
dat.dF = fread(fid, inf, '*single');
dat.dF = reshape(dat.dF, 512,512,[]);
fclose(fid);

%Zscore:
zF = (dat.dF - mean(dat.dF, 3))./std(dat.dF,0,3);
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

%% load data
if ~matches(datanamefluoactivations, datanamefluocurves)
    fid = fopen([DataFolder datanamefluocurves '.dat']);
    dat.dF = fread(fid, inf, '*single');
    dat.dF = reshape(dat.dF, 512,512,[]);
    fclose(fid);
end

fid = fopen([DataFolder datanamehbo '.dat']);
dat.dH = fread(fid, inf, '*single');
dat.dH = reshape(dat.dH,size(dat.dF));
fclose(fid);

fid = fopen([DataFolder datanamehbr '.dat']);
dat.dR = fread(fid, inf, '*single');
dat.dR = reshape(dat.dR,size(dat.dF));
fclose(fid);

clear datanamehbo datanamehbr datanamefluo fid

%% Get mask for brain and movement and outliers
seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'ROImasks_data.mat'], 'img_info');
mask = img_info.logical_mask;

load([DataFolder 'OutlierMask.mat'], 'OutlierFrames', 'OutlierPixels');
OutlierPix = OutlierPixels.HbO + OutlierPixels.HbR + OutlierPixels.hemoCorr_fluo;
OutlierPix(OutlierPix<3) = 0;
OutlierPix(OutlierPix == 3) = 1;

mask = mask.*logical(OutlierPix);
mask = logical(mask);

%take only brain
aF = reshape(aF,[], size(aF,3));
aF = aF(mask(:),:);
dat.dF = reshape(dat.dF,[], size(dat.dF,3));
dat.dF = dat.dF(mask(:),:);
dat.dH = reshape(dat.dH,[], size(dat.dH,3));
dat.dH = dat.dH(mask(:),:);
dat.dR = reshape(dat.dR,[], size(dat.dR,3));
dat.dR = dat.dR(mask(:),:);

%exclude movement and outliers
load([DataFolder 'MovMask.mat'], 'MovMask');
aF = aF.* MovMask;
aF = aF.*OutlierFrames;

Maps.NrOfActs = NaN(512,512);
Maps.NrOfActs(mask) = sum(aF, 2);

clear img_info MovMask seps OutlierMask

datatypes = {'fluo', 'hbo', 'hbr'};
datanames = {'dF', 'dH', 'dR'};

%% get average curves per pixel
for inddata = 1:size(datatypes, 2)
    avcurvesperpixel = NaN(size(aF,1), windowtime);

    for pixel = 1:size(aF, 1)
        % get activations of pixel, but exclude beginning and end (where you
        % cant get a full plot)
        indact = find(aF(pixel,:));
        indact(indact<=(timebefore)) = [];
        indact(indact>=(size(aF,2)-timeafter)) = [];

        if size(indact,2)<1 %if theres no activations
            continue
        end

        temp = NaN(size(indact,2), windowtime);

        for ind2 = 1:size(indact,2)
            activation = indact(ind2); %get curves for all activations of this pixel
            temp(ind2,:) = dat.(datanames{inddata})(pixel,activation-timebefore:activation+timeafter);
        end

        avcurvesperpixel(pixel,:) = mean(temp,1,'omitnan');
    end

    temp = NaN(512*512, windowtime);
    temp(mask,:) = avcurvesperpixel;
    PixelCurves.(datatypes{inddata}) = temp;
    disp([datatypes{inddata} ' done']);
end
clear activation aF fluotemp hbotemp hbrtemp pixel ind2 indact avcurvesperpixel dat datanames inddata temp


%% Maps over brain
Maps.Delay = NaN(512,512);
Maps.Fluopeak = NaN(512, 512);
Maps.HbOpeak = NaN(512, 512);


for indpic = 1:512*512
    if isnan(PixelCurves.fluo(indpic, 1))
        continue
    else
        [Maps.Fluopeak(indpic), fluopeaktiming] = max(PixelCurves.fluo(indpic, :));
        [Maps.HbOpeak(indpic), delay] = max(PixelCurves.hbo(indpic, fluopeaktiming:226));
        Maps.Delay(indpic) = delay-1;
        [Maps.HbRpeak(indpic), ~] = max(PixelCurves.hbr(indpic, fluopeaktiming:226));
    end
end


%% save
save([DataFolder 'NVC_Maps.mat'], "Maps");
clear datanamefluoactivations datanamefluocurves datatypes delay fluopeaktiming imfreq indpic Maps mask OutlierFrames OutlierPix OutlierPixels PixelCurves timeafter timebefore windowtime

%% get rotation of atlas
% load('/media/mbakker/GDrive/P2/GCaMP/BigROI_ReferenceMask.mat', 'AtlasMask');
% atlasmask_ref = AtlasMask;
% 
% seps = strfind(DataFolder, filesep);
% load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask');
% 
% [optimizer, metric] = imregconfig("monomodal");
% movingRegistered = imregister(AtlasMask, atlasmask_ref, "affine", optimizer, metric);



end


