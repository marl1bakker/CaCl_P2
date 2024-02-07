function CorrSpread = SingleSubjectSeedSpread(DataFolder, dataname)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

seps = strfind(DataFolder, filesep);
load([DataFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask', 'regions');

if ~exist([DataFolder dataname], 'file')
    disp([dataname ' not found for ' DataFolder])
    CorrSpread = NaN(size(regions,2), 15);
    return
end

f = waitbar(0, 'Progress SingleSubjectSeedSpread');

% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, inf, '*single');
dat = reshape(dat,512,512,[]);
fclose(fid);
dims = size(dat);
dat = reshape(dat,[],dims(3));

load([DataFolder 'MovMask.mat'], 'MovMask')
dat(dat == 0) = nan;

CorrSpread = NaN(size(regions,2), 15); %15 is a guess, i think we won't have regions that go over 15 loops

for ind = 1:size(regions,2) %go per region
    waitbar(ind/size(regions,2), f)
    if ind == 3 || ind == 8 %if it's the auditory cortex, skip
        continue
    end
    
    Mask = ismember(AtlasMask,ind); %pak nummers van atlas van ind waar je nu bent
    
    %     % get seed and timecourse of it
    %     Tmp = bwmorph(Mask,'shrink',inf); %maak mask 1 pix, noem tmp
    %     Tmp = conv2(Tmp, ones(3),'same')>=1; %maak iets groter
    %     Seed = imerode(Mask, strel('diamond',1)) & Tmp; %krimp de ROI met 1 pixel, pak alleen pixels die ook binnen tmp vallen
    
    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*Mask, [], 1))/sum(Mask(:));
    iY = sum(reshape(Y.*Mask, [], 1))/sum(Mask(:));
    iX = round(iX);
    iY = round(iY);
    
    Seed = zeros(dims(1), dims(2));
    Seed(iY, iX) = 1;
    radius = 3;
    SE = strel('sphere', radius);
    Seed = imdilate(Seed, SE);
    
    %     Orig_seed = Seed; %to have later?
    Seed = logical(Seed);
    timecourse_seed = mean(dat(Seed(:), :),1, 'omitnan');
%     timecourse_seed = dat(iY, iX, :);
    
    periphery = 12;
    per_corrs = [];

    [rho,~] = corr(timecourse_seed', timecourse_seed', 'rows', 'complete'); %rows complete to ignore nan
    per_corr = mean(rho, 'all', 'omitnan');
    per_corrs = [per_corrs; per_corr];
    
    radius = 10;

%     while sum(periphery, 'all') > 11 && radius < 30
    for ind2 = 1:14
        
        %take growing circle around seed
        SE = strel('sphere', radius);
%         radius = radius + 5; % for next round
        periphery = imdilate(Seed, SE);
        periphery = periphery - Seed;
        periphery = periphery.* Mask;
%         imagesc(periphery)
        
        %get timecourse of periphery
        periphery = logical(periphery);
        timecourses_periphery = dat(periphery(:),:);
        
        %take correlation with seed
        [rho, ~] = corr(timecourse_seed', timecourses_periphery', 'rows', 'complete'); 
        per_corr = mean(rho, 'all', 'omitnan');
        per_corrs = [per_corrs; per_corr];
        
        Seed = periphery + Seed;
    end
    
    CorrSpread(ind, 1:size(per_corrs,1)) = per_corrs;
    
end
delete(f)
% figure
% plot(CorrSpread')
end