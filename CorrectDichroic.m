%% Correction for offset red/green channels because of dichroic mirror

function CorrectDichroic(DataFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end
    
%% Settings check-up:
Info = matfile([DataFolder 'AcqInfos.mat']);
AcqInfo = Info.AcqInfoStream;
clear Info;

if( AcqInfo.MultiCam )
    NbIllum = sum(cellfun(@(X) contains(X, 'Illumination'), fieldnames(AcqInfo)));
    Cam1List = {};
    Cam2List = {};
    for ind = 1:NbIllum
        idx = AcqInfo.("Illumination" + int2str(ind)).CamIdx;
        chan = lower(AcqInfo.("Illumination" + int2str(ind)).Color);
        if( contains(chan, 'fluo') )
            chan = ['fluo_' chan(9:11)];    
        end
        if( contains(chan, 'amber') )
            chan = 'yellow';
        end
        if( idx == 1 )
            Cam1List{end+1} = [chan '.dat'];
        else
            Cam2List{end+1} = [chan  '.dat'];
        end 
    end
    clear NbIllum ind idx chan
else
    disp('Only one camera was used. No need to coregister images');
    return;
end

%% Opening Cam1 raw file:
info = matfile([Cam1List{1}(1:(end-3))  'mat']);
fid = fopen(Cam1List{1});
iC1_1 = fread(fid, info.datSize(1,1)*info.datSize(1,2)*10, '*single' );
iC1_1 = reshape(iC1_1, info.datSize(1,1), info.datSize(1,2),[]);
iC1_1 = mean(iC1_1,3);
fid = fopen(Cam1List{2});
iC1_2 = fread(fid, info.datSize(1,1)*info.datSize(1,2)*10, '*single' );
iC1_2 = reshape(iC1_2, info.datSize(1,1), info.datSize(1,2),[]);
iC1_2 = mean(iC1_2,3);

if( mean(iC1_1(:)) > 1e3 )
    iC1_1 = iC1_1./mean(iC1_1(:));
else
    iC1_1 = zeros(size(iC1_1));
end
if( mean(iC1_2) > 1e3 )
    iC1_2 = iC1_2./mean(iC1_2(:));
else
    iC1_2 = zeros(size(iC1_1));
end
iC1 = iC1_1 + iC1_2;
if(mean(iC1(:)) < 0.1 )
    disp('Signals recorded on camera 1 are too low to perform coregistration.')
    return;
end
iC1 = iC1./mean(iC1(:));
clear iC1_* 

%% Opening Cam2 raw file:
fid = fopen(Cam2List{1});
iC2_1 = fread(fid, info.datSize(1,1)*info.datSize(1,2)*10, '*single' );
iC2_1 = reshape(iC2_1, info.datSize(1,1), info.datSize(1,2),[]);
iC2_1 = mean(iC2_1,3);
fid = fopen(Cam2List{2});
iC2_2 = fread(fid, info.datSize(1,1)*info.datSize(1,2)*10, '*single' );
iC2_2 = reshape(iC2_2, info.datSize(1,1), info.datSize(1,2),[]);
iC2_2 = mean(iC2_2,3);

if( mean(iC2_1(:)) > 1e3 )
    iC2_1 = iC2_1./mean(iC2_1(:));
else
    iC2_1 = zeros(size(iC2_1));
end
if( mean(iC2_2(:)) > 1e3 )
    iC2_2 = iC2_2./mean(iC2_2(:));
else
    iC2_2 = zeros(size(iC2_1));
end
iC2 = iC2_1 + iC2_2;
if(mean(iC2(:)) < 0.1 )
    disp('Signals recorded on camera 2 are too low to perform coregistration.')
    return;
end
% iC2 = iC2./mean(iC2(:));
%just take the red, not the f'ed up fluo. Increase the small differences
iC2_1 = (iC2_1 - min(iC2_1(:)))./(max(iC2_1(:))-min(iC2_1(:)));
iC2 = adapthisteq(iC2_1);

clear iC2_* fid 
%% Coregistration parameters:

% ****WARNING**** 
% These following lines are only in the specific case where camera 2 image
% is flipped (mirrored) compare to camera 1 image.
% To be executed only in that case!!!
% InitialT =affine2d( [-1 0 0; 0 1 0; size(iC2,1) 0 1]);
%Value of InitialT when there is no flip:
InitialT = affine2d([1 0 0; 0 1 0; 0 0 1]);
%end of WARNING

[opt, met] = imregconfig("multimodal");
opt.GrowthFactor = 1.05;
% opt.GrowthFactor = 2;

opt.Epsilon = 1.5e-6;

% opt.InitialRadius = 2.5e-3; %old
opt.InitialRadius = 1e-4;

% opt.MaximumIterations = 100; %old
opt.MaximumIterations = 200;

%% Coregistration 
figure(1); 
imshowpair(imwarp(iC2, InitialT, 'OutputView', imref2d(size(iC2))), iC1);
tform = imregtform(imwarp(iC2, InitialT, 'OutputView', imref2d(size(iC2))), iC1, 'similarity', opt, met);
tform.T = InitialT.T*tform.T;

figure(2);
imshowpair(imwarp(iC2,tform,'OutputView', imref2d(size(iC1))), iC1);

figure(3)
for ind = 1:30
    imagesc(iC1)
    pause(0.2)
    imagesc(imwarp(iC2,tform,'OutputView', imref2d(size(iC1))))
    pause(0.2)
end

% figure(3)
% for ind = 1:50
%     imagesc(iC1)
%     pause(0.2)
%     imagesc(iC2)
%     pause(0.2)
% end

answer = questdlg('Was the coreg good?',...
    'Coregistration cameras',...
    'Yes', 'No', 'No');

switch answer
    case 'Yes'
        % :)
        disp('overwriting .dat files')
        close figure 1
        close figure 2
        close figure 3
    case 'No'
        disp(['did not coregister recording ' DataFolder])
        close figure 1
        close figure 2
        close figure 3
        return;
end

%% Correction:
%If the last image shown was coregistered properly on the camera 1 image,
%execute this part to correct camera 2 images:

for ind = 1:size(Cam2List,2)
    fid = fopen(Cam2List{ind});
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, info.datSize(1,1), info.datSize(1,2),[]);
    dat = imwarp(dat, tform, 'OutputView', imref2d(info.datSize));

    fclose(fid);
    fid = fopen(Cam2List{ind},'w');
    fwrite(fid,dat,'single');
    fclose(fid);
end

%% If you want to keep the transformation matrix :
% save('tform.mat', 'tform');
    
end






























%% old
% if( ~strcmp(DataFolder(end), filesep) )
%     DataFolder = [DataFolder filesep];
% end
%     
% %load single picture data
% fid = fopen([DataFolder 'green.dat']);
% green = fread(fid, 512*512, '*single');
% green = reshape(green, 512, 512);
% fclose(fid);
% fid = fopen([DataFolder 'red.dat']);
% red = fread(fid, 512*512, '*single');
% red = reshape(red, 512,512);
% fclose(fid);
% % imshowpair(green, red)
% 
% % To increase small characteristics
% green = (green - min(green(:)))./(max(green(:))-min(green(:)));
% red = (red - min(red(:)))./(max(red(:))-min(red(:)));
% red = adapthisteq(red);
% green = adapthisteq(green);
% 
% % imshowpair(green, red)
% [opt, met] = imregconfig('multimodal');
% 
% tform = imregtform(green, red,'similarity',opt, met); %moving pic is green
% % imshowpair(imwarp(green,tform, 'OutputView', imref2d(size(green))), red)
% 
% 
% % Load all data to transform on -- green
% fid = fopen([DataFolder 'green.dat']);
% green = fread(fid, inf, '*single');
% green = reshape(green,512,512,[]);
% fclose(fid);
% 
% for ind = 1:size(green,3)
%     green(:,:,ind) = imwarp(green(:,:,ind), tform, 'OutputView', imref2d([512 512]));
% end
%     
% % save([DataFolder 'green.mat'], green);
% fid = fopen([DataFolder 'green.dat'],'w');
% fwrite(fid, green,'single');
% fclose(fid);
% 
% % Load all data to transform on -- fluo
% fid = fopen([DataFolder 'fluo_567.dat']);
% fluo = fread(fid, inf, '*single');
% fluo = reshape(fluo,512,512,[]);
% fclose(fid);
% 
% for ind = 1:size(fluo,3)
%     fluo(:,:,ind) = imwarp(fluo(:,:,ind), tform, 'OutputView', imref2d([512 512]));
% end
%     
% % save([DataFolder 'green.mat'], green);
% fid = fopen([DataFolder 'fluo.dat'],'w');
% fwrite(fid, fluo,'single');
% fclose(fid);
% 
% end
% 
% 


