%P2
% This function is used to coregister different acquisitions of the same
% mouse with each other, so that the brain is at the same place. It needs
% the function CoregistrationManual. 

% path_fixed is the folder of the mouse with data to base your fixation on,
% in this case A1. path_moving are the datafolders with the files that you
% want to coregister, in this case A2 and A3. dataname is the file you base
% your coregistration on, in this case green or speckle. ToBeCorrected are
% the files that you want to move. In case of GCaMP, that is green, red,
% and fluo, in case of speckle it is speckle. Give like {'green', 'red', 'fluo_567.dat'}

function Coregistration_speckle(path_fixed, path_moving, dataname, ToBeCorrected)

if ~exist('dataname', 'var')
    dataname = 'BFI.dat';
elseif size(dataname,2)<4 || ~matches(dataname(end-3:end), '.dat')
    dataname = [dataname '.dat'];
end

if( ~strcmp(path_moving(end), filesep) )
    path_moving = [path_moving filesep];
end
if( ~strcmp(path_fixed(end), filesep) )
    path_fixed = [path_fixed filesep];
end

%% Step 1: Load data
load([path_fixed 'meanBFI.mat'], 'meanBFI');
datfixed = meanBFI;
datfixed(datfixed>1500) = 1500;

load([path_moving 'meanBFI.mat'], 'meanBFI');
datmov = meanBFI;
datmov(datmov>1500) = 1500;

% load([path_fixed dataname(1:end-4) '.mat'], 'datSize');
% 
% fid = fopen([path_fixed dataname]);
% dat = fread(fid, datSize(1)*datSize(2), '*single');
% dat = reshape(dat, datSize(1),datSize(2));
% fclose(fid);
% datfixed = dat;
% 
% fid = fopen([path_moving dataname]);
% dat = fread(fid, datSize(1)*datSize(2), '*single');
% dat = reshape(dat, datSize(1),datSize(2));
% fclose(fid);
% datmov = dat;
% clear dat

%% Step 3: Prep, make differences more clear
% datfixed = datfixed./mean(datfixed(:));
% datfixed = (datfixed-min(datfixed(:)))./(max(datfixed(:))-min(datfixed(:)));
% datfixed(datfixed < 0) = 0;
% datfixed(datfixed > 1) = 1;
% datfixed = adapthisteq(datfixed);
% 
% datmov = datmov./mean(datmov(:));
% datmov = (datmov-min(datmov(:)))./(max(datmov(:))-min(datmov(:)));
% datmov(datmov < 0) = 0;
% datmov(datmov > 1) = 1;
% datmov = adapthisteq(datmov);

%% Step 4: CoReg
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 1e6;
optimizer.MinimumStepLength = 1e-5;
optimizer.MaximumStepLength = 1e-2;
optimizer.GradientMagnitudeTolerance = 1e-3;
tform = imregtform(datmov, datfixed, 'affine', optimizer, metric,...
    'DisplayOptimization', false, 'PyramidLevels', 3);

%% Step 5:Confirmation
f1 = figure;
f1.Name = 'old';
imshowpair(datfixed, datmov);
dat_corr = imwarp(datmov,tform,'OutputView',imref2d(size(datfixed)));
f2 = figure;
f2.Name = 'new';
imshowpair(datfixed, dat_corr);
% saveas(gcf,[path_moving 'CoregistrationComp.png']);
% 
pause 
answer = questdlg('Does it make sense?', ...
	'Coregistration', ...
	'Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        disp('≧◠‿●‿◠≦    ᕙ(^▿^-ᕙ)');
    case 'No'
        disp('Coreg. not applied')
%         return
%         close all % BAD PROGRAMMING CHANGE LATER
        close(f1, f2)
        error('Coregistration not good enough to apply')
end
% close all % BAD PROGRAMMING CHANGE LATER
close(f1, f2)
disp('applying coregistration...')

%% Step 5 B: Save/load
% save([path_moving 'tform.mat'],'tform');

%% Step 6: Apply to all images
% if( bApply )
if ~exist('ToBeCorrected', 'var')
    ToBeCorrected = {'green.dat', 'red.dat', 'fluo_567.dat'};
end

for index = 1:size(ToBeCorrected,2)
    dataname = ToBeCorrected{index};
    fid = fopen([path_moving dataname]);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, datSize(1),datSize(2),[]);
    fclose(fid);
    
    for ind = 1:size(dat,3)
        dat(:,:,ind) = imwarp(squeeze(dat(:,:,ind)),tform,'OutputView',imref2d(size(dat(:,:,1))),'interp','nearest');
    end
    
    fid = fopen([path_moving dataname], 'w');
    fwrite(fid,dat,'*single');
	fclose(fid);   
end


end %of function
