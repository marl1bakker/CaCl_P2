%P2
% path_fixed =  '/media/mbakker/GDrive/P2/GCaMP/M13/A1-R2/CtxImg';
% path_moving = '/media/mbakker/GDrive/P2/GCaMP/M13/A2-R1/CtxImg';
function Coregistration(path_fixed, path_moving, dataname)

if ~exist('dataname', 'var')
    dataname = 'green.dat';
end

if( ~strcmp(path_moving(end), filesep) )
    path_moving = [path_moving filesep];
end
if( ~strcmp(path_fixed(end), filesep) )
    path_fixed = [path_fixed filesep];
end

%% Step 1: Load data
%check if both acquisitions have done coregistration within, otherwise the
%orientation is gonna be different
% if( exist([path_fixed 'IntraCoReg.mat'], 'file') ) && ( exist([path_moving 'IntraCoReg.mat'], 'file') )
%     disp('both paths have been coregistered within their own acquisition')
% else
%     disp('one of the acquisitions have not been coregistered within the acquisition')
%     datH12 = 0;
%     return 
% end

fid = fopen([path_fixed dataname]);
dat = fread(fid, 512*512, '*single');
dat = reshape(dat, 512,512);
fclose(fid);
datfixed = dat;

fid = fopen([path_moving dataname]);
dat = fread(fid, 512*512, '*single');
dat = reshape(dat, 512,512);
fclose(fid);
datmov = dat;
clear dat

%% Step 3: Prep, make differences more clear
datfixed = datfixed./mean(datfixed(:));
datfixed = (datfixed-min(datfixed(:)))./(max(datfixed(:))-min(datfixed(:)));
datfixed(datfixed < 0) = 0;
datfixed(datfixed > 1) = 1;
datfixed = adapthisteq(datfixed);

datmov = datmov./mean(datmov(:));
datmov = (datmov-min(datmov(:)))./(max(datmov(:))-min(datmov(:)));
datmov(datmov < 0) = 0;
datmov(datmov > 1) = 1;
datmov = adapthisteq(datmov);

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
        close all % BAD PROGRAMMING CHANGE LATER
        error('Coregistration not good enough to apply')
end
close all % BAD PROGRAMMING CHANGE LATER
disp('applying coregistration...')
%% Step 5 B: Save/load
% save([path_moving 'tform.mat'],'tform');

%% Step 6: Apply to all images
% if( bApply )
ToBeCorrected = {'green.dat', 'red.dat', 'fluo_567.dat'};

for index = 1:size(ToBeCorrected,2)
    dataname = ToBeCorrected{index};
    fid = fopen([path_moving dataname]);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, 512,512,[]);
    fclose(fid);
    
    for ind = 1:size(dat,3)
        dat(:,:,ind) = imwarp(squeeze(dat(:,:,ind)),tform,'OutputView',imref2d(size(dat(:,:,1))),'interp','nearest');
    end
    
    fid = fopen([path_moving dataname], 'w');
    fwrite(fid,dat,'*single');
	fclose(fid);   
end
% else
%     
% end

end %of function
