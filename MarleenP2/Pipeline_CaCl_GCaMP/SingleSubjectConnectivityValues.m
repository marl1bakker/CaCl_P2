function [ConnValues, nrpixels] = SingleSubjectConnectivityValues(DataFolder, dataname, threshold)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

%check if you have the clustered ROI, load it
idx = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:idx(end-2));

if ~exist([pathROI 'BigROI.mat'], 'file')
    disp(['BigROI does not exist for ' DataFolder])
    return
end

if ~exist('threshold', 'var')
    threshold = 0.6;
end

% Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
% Acq = DataFolder(idx(end-2)+1:idx(end-1)-1);

clear idx
f = waitbar(0, 'Progress SingleSubjectConnectivityValues');

%% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, inf, '*single');
fclose(fid);
dat = reshape(dat, 512*512, []);

load([DataFolder 'MovMask.mat'], 'MovMask')
dat = dat .* MovMask; 
dat(dat == 0) = NaN;

load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');
GenMask = logical(AtlasMask);
datbrain = dat(GenMask(:),:);
dat = reshape(dat, 512,512, []);
dims = size(dat);

ConnValues = nan(1,size(regions,2));
nrpixels = sum(~isnan(datbrain(:,find(MovMask, 1))), 'all');

%% Do everything per seedname you gave in
for ind = 1:size(regions, 2)
    waitbar(ind/size(regions,2), f)
%     disp(regions{ind})
    
    %% Get middle of ROIs
    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iY = sum(reshape(Y.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iX = round(iX);
    iY = round(iY);
    
    if isnan(iX) || isnan(iY)
        continue
    end
    
    %% Calculate the seed pixel correlation map
    Seeddat = dat(iY, iX, :);
    Seeddat = reshape(Seeddat, 1, []);
    
    % without move
    [rho, ~] = corr(Seeddat', datbrain(:,:)', 'rows', 'pairwise');

    ConnValues(ind) = sum(rho>threshold, 'all');
    
end

delete(f)
end