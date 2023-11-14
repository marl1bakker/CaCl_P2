%test mice
%% gcamp
% M1 check for imagesclassification if it works
cd '/media/mbakker/SSD-2TB/Marleen/M1-W1-BR2'

%datafolder, savefolder, binning (with BR2 it's 2x2, otherwise its 1x1), 
ImagesClassification(pwd, pwd, 1, 1, 0, 0, 1);


%% speckle
cd '/media/mbakker/PJM - HDD - 2/Marleen/Speckle/M1-A2-SR1'
ImagesClassification(pwd, pwd, 1, 1, 1, 0, 'Internal-main');
Ana_Speckle(pwd,0)

%% test data transfer
cd '/media/mbakker/PJM - HDD - 2/Marleen/transfertest'
ImagesClassification(pwd, pwd, 1, 1, 1, 0, 'Internal-main');

cd '/media/mbakker/PJM - HDD - 2/Marleen/transfertestnotransfer'
ImagesClassification(pwd, pwd, 1, 1, 1, 0, 'Internal-main');

% transfering data simultaneously does not make a huge difference

%% Check if acquisitions are going well. Take acquisition 3 of M7, GCaMP
Pipeline_CaCl_GCaMP('/media/mbakker/data1/P2/M7-A3-R1', 1)

% fid = fopen('green.dat');
% datgreen = fread(fid, inf, '*single');
% datgreen = reshape(datgreen, 512,512,[]);

datgreen = NormalisationFiltering(pwd, 'green', 1/120, 0.1, 1,0);

for ind = 1:1000
    imagesc(datgreen(:,:,ind), [0.99,1.01]);
    title(num2str(ind));
    pause(0.0025);
end

% fid = fopen('red.dat');
% datred = fread(fid, inf, '*single');
% datred = reshape(datred, 512,512,[]);

datred = NormalisationFiltering(pwd, 'red', 1/120, 0.1, 1,0);

for ind = 1:1000
    imagesc(datred(:,:,ind), [0.99, 1.01]);
    title(num2str(ind));
    pause(0.0025);
end


fid = fopen('fluo_567.dat');
datfluo = fread(fid, inf, '*single');
datfluo = reshape(datfluo, 512,512,[]);

imagesc(datfluo(:,:,2));

for ind = 1:500
    imagesc(datfluo(:,:,ind), [0.2,20000]);
    title(num2str(ind));
    pause(0.05);
end

% OutData = NormalisationFiltering(FolderData, FileData, lowFreq,...
%     highFreq, bDivide, bExpfit, varargin)
datfluo = NormalisationFiltering(pwd, 'fluo_567', 0.3, 3, 1,0);

for ind = 1:500
    imagesc(datfluo(:,:,ind), [0.95,1.05]);
    title(num2str(ind));
    pause(0.05);
end


%%
cd '/media/mbakker/data1/P2/M7-A3-RS1'
% ImagesClassification(pwd, pwd, 1, 1, 1, 0, 'Internal-main');

datspeckle = NormalisationFiltering(pwd, 'speckle', 1/120, 5, 1,0);

kernel = ones(1,1,5,'single'); 

pix2 = 0;
for ind = 0:7
    pix1 = pix2 + 1;
    pix2 = pix1+127;
    datspeckle(pix1:pix2,:,:) = stdfilt(datspeckle(pix1:pix2,:,:), kernel);
    disp(num2str(ind));
end

for ind = 1:200
    imagesc(-log10(datspeckle(:,:,ind)),[1.5,2.5]);
    title(num2str(ind));
    pause(0.05);
end

%% check running wheel
aiFilesList = dir([pwd filesep 'ai_*.bin']);

% Opening of the files:
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([pwd filesep aiFilesList(ind).name],...
        'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, 11, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],11);
    AnalogIN = [AnalogIN; tmp];
end

plot(AnalogIN(:,5));

%%
fid = fopen('green.dat');
datgreen = fread(fid, inf, '*single');
datgreen = reshape(datgreen, 512,512,[]);

%% image class 
[Mice, AcqList] = MakeAcqList;

for ind =1:size(AcqList,2)
    Pipeline_CaCl_GCaMP(AcqList(ind).name,1)
end