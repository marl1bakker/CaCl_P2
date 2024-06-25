%% Mark frames with movement
% 1 is no movement, 0 is movement

function MarkMovedFrames(DataFolder, SaveFolder)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('SaveFolder', 'var')
    SaveFolder = DataFolder;
end

if( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = [SaveFolder filesep];
end

if exist([SaveFolder 'MovMask.mat'], 'file')
    return
end

if contains(SaveFolder, 'Speckle')
    Info = open([SaveFolder 'speckle.mat']);
elseif contains(SaveFolder, 'GCaMP')    
    Info = open([SaveFolder 'fluo_567.mat']);
else
    disp('Type of aquisition not recognized!')
%     return
end

acqtime = floor(Info.datLength/Info.Freq); %sec
nrofframes = Info.Freq*acqtime;

% Opening of the files:
aiFilesList = dir([DataFolder 'ai_*.bin']);
AnalogIN = [];
for ind = 1:size(aiFilesList,1)
    data = memmapfile([DataFolder aiFilesList(ind).name],...
        'Offset', 5*4, 'Format', 'double', 'repeat', inf);
    tmp = data.Data;
    tmp = reshape(tmp, 1e4, 11, []);
    tmp = permute(tmp,[1 3 2]);
    tmp = reshape(tmp,[],11);
    AnalogIN = [AnalogIN; tmp];
end

TreadmillFull = AnalogIN(:,5);
% TreadmillFull = TreadmillFull(1:6000000);
%the analog in starts recording before the actual acquisition. We need AI1
%to see where the acquisition started and how much of the AI channel is
%from before the acquisition

CamTrigger = AnalogIN(:,1);
CamTrigger(CamTrigger > 4) = 5;
CamTrigger(CamTrigger < 1) = 0;
if CamTrigger(1) == 5
    firstzero = find(CamTrigger == 0, 1);
    StartAcq = firstzero + find(CamTrigger(firstzero:end) == 5, 1);
else 
    StartAcq = find(CamTrigger == 5, 1);
end

TreadmillFull = TreadmillFull(StartAcq:end);

clear aiFilesList AnalogIN data ind tmp firstzero 

% center it around 0 so you can get the maximum outlier easier (if not,
% your max outlier might be below the standard of about 1.64, and it
% will be difficult to detect it)
TreadmillFull = TreadmillFull - median(TreadmillFull);

%get one treadmill value per frame, make values absolute
Treadmill = [];
ind1 = 1;
for ind = 1:nrofframes
%     ind2 = round(6000000/nrofframes * ind); %i think this is based on acqtime, but did not check
    ind2 = round((acqtime*10000)/nrofframes * ind);
    Treadmill(end+1) = max(abs(TreadmillFull(ind1:ind2)), [], 1);
    ind1 = ind2+1;
end

clear TreadmillFull ind ind1 ind2

%Tried to do it based on outliers, but think cutoff is better
%     Outliers = isoutlier(Treadmill, 'gesd'); %checked, seems like gesd gives best results
cutoff = 0.01;
temp = Treadmill < cutoff;

%now the problem is that it is missing some values in the middle, where
%there was mov but the value was around the median.
% expand the movement frames by 10. Not the most elegant way but it works.
%it now marks things with no movement as 1 and things with mov as 0
MovMask = ones(size(temp)); % if you want 0 for no mov en 1 for mov, change this to zeros(size(temp)) and ...
% buffer = round(Info.Freq/3)*2;
buffer = Info.Freq; % 1 sec buffer
for ind = 1:size(temp,2)
    if temp(ind) == 0
        if ind-buffer <= 0 
            ind1 = 1;
            ind2 = ind+buffer;
        elseif ind+buffer > size(temp,2)
            ind1 = ind-buffer;
            ind2 = size(temp,2);
        else
            ind1 = ind-buffer;
            ind2 = ind+buffer;
        end
        MovMask(ind1:ind2) = 0; % cont. ... this part to MovMask(ind1:ind2) = 1;
    end
end

% if Info.datLength > 9000
if Info.datLength > (Info.Freq*acqtime)
%     nrextraframes = Info.datLength - 9000;
    nrextraframes = Info.datLength - (Info.Freq*acqtime);
    MovMask(end+nrextraframes) = 0; %make extra frames 0 so you mark them as movement
end

clear temp ind ind1 ind2 cutoff

% figure()
% plot(MovMask)
% hold on
% plot(Treadmill)
% pause()

close all

save([SaveFolder 'MovMask.mat'], 'MovMask');

end