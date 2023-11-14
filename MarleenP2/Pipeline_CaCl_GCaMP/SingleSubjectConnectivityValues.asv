function [ConnValues, nrpixels] = SingleSubjectConnectivityValues(DataFolder, dataname)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

%check if you have the clustered ROI, load it
idx = strfind(DataFolder, filesep); %zoek alle plekken van fileseps in naam
pathROI = DataFolder(1:idx(end-2));

if ~exist([pathROI 'BigROI.mat'], 'file')
    disp(['BigROI does not exist for ' DataFolder])
    return
end

% Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
% Acq = DataFolder(idx(end-2)+1:idx(end-1)-1);

clear idx

%% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512,512, []);

load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');
GenMask = logical(AtlasMask);
dat = dat.* GenMask;
dat(dat == 0) = NaN;
nrpixels = sum(~isnan(dat(:,:,1)),'all');

dims = size(dat);
ConnValues = nan(1,size(regions,2));

sum(~isnan(dat(:,:,1)), 'all'); %how many pixels do you have for this mousebrain
%% Do everything per seedname you gave in
for ind = 1:size(regions, 2)
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
    dat = reshape(dat,dims);
    Seeddat = dat(iY, iX, :);
    Seeddat = reshape(Seeddat, 1, []);
    dat = reshape(dat, dims(1)*dims(2), []);
    [rho, ~] = corr(Seeddat', dat(:,:)');
    rho = reshape(rho, dims(1), dims(2));
    ConnValues(ind) = sum(rho>0.6, 'all');
    
end
end