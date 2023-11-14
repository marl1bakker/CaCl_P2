function CoregistrationManual(path_fixed, path_moving, dataname)

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


%% Do manual Coregistration
tform = AlignFrame_Manual(datfixed, datmov);

%% Apply to all images
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
    disp([dataname ' saved'])
end

end %of function
