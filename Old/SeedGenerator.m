function SeedGenerator(DataFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

% load ROI
fileseps = strfind(DataFolder, filesep);
load([DataFolder(1:fileseps(end-2)) 'ROImasks_data.mat']);

for ind = 1:size(ROI_info,2)
    ROI = ROI_info(ind).Stats.ROI_binary_mask;
    %% Find Centroid
    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:512, 1:512);
    iX = sum(X(:).*ROI(:))/sum(ROI(:));
    iY = sum(Y(:).*ROI(:))/sum(ROI(:));
    iX = round(iX);
    iY = round(iY);
    singlepixel = [iX iY];
    %expand slightly
    Seed = zeros(512);
    Seed(iY, iX) = 1;
    Seed = conv2(Seed, fspecial('disk',3)>0,'same')>0;
    
    [ROI_info(ind).('Seeds')] = Seed;
    [ROI_info(ind).('PixelSeeds')] = singlepixel;
end

save([DataFolder(1:fileseps(end-2)) 'ROImasks_data.mat'], 'ROI_info', 'img_info');
end