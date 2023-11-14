% choose one pixel in the brain, and make a correlationmap of that. 
% ManualInput 0 - dont ignore existing files
% ManualInput 1 - overwrite existing files

%% P2
% Input Seedname can be several, give like {'seed1', 'seed2'} etc.
%{'VisualROI_R'} SensoryROI_R AuditoryROI_R UnknownROI_R MotorROI_R

%goes with BigROI standard, because the small rois will give way too many
%points. 

function Allrho = SingleSubjectSPCM(DataFolder, dataname, Overwrite, GSR)

if ~exist('Overwrite', 'var')
    Overwrite = 0;
end

if ~exist('GSR', 'var')
    GSR = 0;
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

idx = strfind(DataFolder, filesep); 
pathROI = DataFolder(1:idx(end-2));
if ~exist([pathROI 'BigROI.mat'], 'file')
    disp(['BigROI does not exist for ' DataFolder])
    return
end

Mouse = DataFolder(idx(end-3)+1:idx(end-2)-1);
Acq = DataFolder(idx(end-2)+1:idx(end-1)-1);

if GSR == 0 && Overwrite == 0 && exist(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM.mat'], 'file')
    disp(['SPCM ' Mouse ' ' Acq 'already made'])
    return
elseif GSR == 1 && Overwrite == 0 && exist(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM-GSR.mat'], 'file')
    disp(['SPCM ' Mouse ' ' Acq 'already made'])
    return
end

clear idx 

%% Get data
fid = fopen([DataFolder dataname]);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512,512, []);

load([pathROI 'BigROI.mat'], 'AtlasMask', 'BigROI', 'regions');
GenMask = logical(AtlasMask);
dat = dat.* GenMask;
dat(dat == 0) = NaN;

dims = size(dat);
Allrho = nan(size(regions, 2), 512, 512);

%% GSR
if GSR == 1
    dat = reshape(dat,[], dims(3));
    mS = mean(dat,1, 'omitnan');
    
    X = [ones(size(mS)); mS];
    B = X'\dat';
    A = (X'*B)';
    dat = dat./A;
    dat = reshape(dat,dims);
    clear h Mask mS X B A;
end

%% make the SPCM
    % Do everything per seedname you gave in
for ind = 1:size(regions, 2) 
    disp(regions{ind})
    
    % Get middle of ROIs
    % Get centroid of ROI based on weight
    [X, Y] = meshgrid(1:dims(1), 1:dims(2));
    iX = sum(reshape(X.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iY = sum(reshape(Y.*BigROI.(regions{ind}), [], 1))/sum(BigROI.(regions{ind})(:));
    iX = round(iX);
    iY = round(iY);
    
    % Calculate the seed pixel correlation map
    dat = reshape(dat,dims);
    Seeddat = dat(iY, iX, :);
    Seeddat = reshape(Seeddat, 1, []);
    dat = reshape(dat, dims(1)*dims(2), []);
    [rho, ~] = corr(Seeddat', dat(:,:)');
    rho = reshape(rho, dims(1), dims(2));
    
    Allrho(ind, :,:) = rho; 
end

if GSR == 0
    save(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM.mat'], 'Allrho')
elseif GSR == 1
    save(['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq '-SPCM-GSR.mat'], 'Allrho')
end

%     %% plot
%     if plot == 1
%         figure()
%         if GSR == 0
%             imagesc(rho, [0 1]) %without gsr
%         else
%             imagesc(rho, [-1 1])
%         end
%         %     colormap jet
%         load('/media/mbakker/data1/Hypoxia/SeedPixelCorrMap/NL.mat');
%         colormap(NL)
%         colorbar
%         
%         axis image
%         
%         line([5 55], [5 5], 'color', 'yellow'); %50 pix per mm
%         line([5 5], [5 55], 'color', 'yellow');
%         
%         if GSR == 1
%             saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
%                 '_' regions{ind} '_GSR.tiff'], 'tiff');
%             saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
%                 '_' regions{ind} '_GSR.eps'], 'epsc');
%         else
%             saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
%                 '_' regions{ind} '_noGSR.tiff'], 'tiff');
%             saveas(gcf, ['/media/mbakker/GDrive/P2/GCaMP/SPCM/Individual/' Mouse '-' Acq ...
%                 '_' regions{ind} '_noGSR.eps'], 'epsc');
%         end
%         close gcf
%     end
    
end