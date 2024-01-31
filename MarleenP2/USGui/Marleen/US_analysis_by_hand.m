%% Analysis by hand, US
load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS.mat');

%% M33-R
% PW:
pwv_file='C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\Temp\2023-03-31-15-30-26_M33 PW mode 1 right-2023-03-31-13-12-26_1.avi';
[mean_velocity, PI, RI] = US_analysis_by_hand_PW(pwv_file);
% %Results were: 
% mean_velocity = 236.6279;
% PI = 2.0957;
% RI = 0.7990;

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M33-R"));
ResultsUS.MeanVelocity(mouseindex) = mean_velocity;
ResultsUS.PulsatilityIndex(mouseindex) = PI;
ResultsUS.ResistenceIndex(mouseindex) = RI;
clear pwv_file mean_velocity PI RI mouseindex
%% M32-L
% MMode:
mmode_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M32\xx - 2023-03-31-15-30-26_M32 M mode 3 left-2023-03-31-11-11-06_1.avi";
[diameter_change, pct_diameter_change] = US_analysis_by_hand_MMode(mmode_file);

% diameter_change = 0.10985;
% pct_diameter_change = 0.25727;

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M32-L"));
ResultsUS.DiameterChange(mouseindex) = diameter_change;
ResultsUS.PCTDiameterChange(mouseindex) = pct_diameter_change;

%% M32-L
% PWMode:
pwv_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M32\xx - 2023-03-31-15-30-26_M32 PW mode 6 left-2023-03-31-11-24-44_1.avi";
[mean_velocity, PI, RI] = US_analysis_by_hand_PW(pwv_file);
% mean_velocity = 148.8752;
% PI = 2.7281;
% RI = 0.84251;

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M32-L"));
ResultsUS.MeanVelocity(mouseindex) = mean_velocity;
ResultsUS.PulsatilityIndex(mouseindex) = PI;
ResultsUS.ResistenceIndex(mouseindex) = RI;
clear pwv_file mean_velocity PI RI mouseindex
%% M31-L 
% MMode:
mmode_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M31\xxx - 2023-03-31-15-30-26_M31 M mode 1 left-2023-03-31-10-00-05_1.avi";
[diameter_change, pct_diameter_change] = US_analysis_by_hand_MMode(mmode_file);
% diameter_change = 0.072;
% pct_diameter_change = 0.22773;
%note: difficult to get 2 heartbeats uninterrupted

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M31-L"));
ResultsUS.DiameterChange(mouseindex) = diameter_change;
ResultsUS.PCTDiameterChange(mouseindex) = pct_diameter_change;
%% M29-R 
% MMode
mmode_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M29\Used - 2023-03-12-15-00-18_M29 M Mode 3 12-3-2023-03-12-13-47-00_1.avi";
[diameter_change, pct_diameter_change] = US_analysis_by_hand_MMode(mmode_file);

% diameter_change = 0.081692;
% pct_diameter_change = 0.1866;
% note: cutoffs changed

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M29-R"));
ResultsUS.DiameterChange(mouseindex) = diameter_change;
ResultsUS.PCTDiameterChange(mouseindex) = pct_diameter_change;
%% M26-R
% MMode
mmode_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M26\Used - 2023-03-12-15-00-18_M26 M Mode 3 12-3-2023-03-12-14-50-45_1.avi";
[diameter_change, pct_diameter_change] = US_analysis_by_hand_MMode(mmode_file);

% diameter_change = 0.081231;
% pct_diameter_change = 0.12213;

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M26-R"));
ResultsUS.DiameterChange(mouseindex) = diameter_change;
ResultsUS.PCTDiameterChange(mouseindex) = pct_diameter_change;
clear mmode_file diameter_change pct_diameter_change mouseindex
%% M25-R
% PW mode
pwv_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M25\Used - 2023-03-12-15-00-18_M25 PW Mode 2 12-3-2023-03-12-11-18-41_1.avi";
[mean_velocity, PI, RI] = US_analysis_by_hand_PW(pwv_file);

% mean_velocity = 155.1573;
% PI = 3.0344;
% RI = 0.85527;

% get mouse:
mouseindex = find(strcmp(ResultsUS.Row, "M25-R"));
ResultsUS.MeanVelocity(mouseindex) = mean_velocity;
ResultsUS.PulsatilityIndex(mouseindex) = PI;
ResultsUS.ResistenceIndex(mouseindex) = RI;
clear pwv_file mean_velocity PI RI mouseindex

%% Dont forget to save!
save('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\ResultsUS+byhand.mat', 'ResultsUS');


















%% Analysis by hand, US
% pwv_file='C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\Marleen\Temp\2023-03-31-15-30-26_M33 PW mode 1 right-2023-03-31-13-12-26_1.avi';

function [mean_velocity, PI, RI] = US_analysis_by_hand_PW(pwv_file)
%% set up
% Adjusted from displayImages function
pwv_v=VideoReader(pwv_file);
% 3 frames
v.CurrentTime=1-1/30;
tmp=readFrame(pwv_v);
app.pwv_frames = []; %to make sure the pwv_frames is in 4D double and not uint
app.pwv_frames(1:size(tmp,1),:,:,1)=tmp;
% app.PWVFrameIndexDropDown.Items = [app.PWVFrameIndexDropDown.Items num2str(1)];
app.current_results.pwv_mode(1).use_data=0;

if pwv_v.Duration >1 && pwv_v.Duration <2
    pwv_v.CurrentTime=1+1/30;
    tmp=readFrame(pwv_v);
    app.pwv_frames(1:size(tmp,1),:,:,2)=tmp;
    % app.PWVFrameIndexDropDown.Items = [app.PWVFrameIndexDropDown.Items num2str(2)];
    app.current_results.pwv_mode(2).use_data=0;
end

if pwv_v.Duration>=2
    % if v.Duration>=1+20/30
    pwv_v.CurrentTime=2;
    tmp=readFrame(pwv_v);
    app.pwv_frames(1:size(tmp,1),:,:,2)=tmp;
    % app.PWVFrameIndexDropDown.Items = [app.PWVFrameIndexDropDown.Items num2str(2)];
    app.current_results.pwv_mode(2).use_data=0;
end

if pwv_v.Duration>=3+1/30
    pwv_v.CurrentTime=3+1/30;
    tmp=readFrame(pwv_v);
    app.pwv_frames(1:size(tmp,1),:,:,3)=tmp;
    % app.PWVFrameIndexDropDown.Items = [app.PWVFrameIndexDropDown.Items num2str(3)];
    app.current_results.pwv_mode(3).use_data=0;
end

%% start
%this while loop is to make sure that if the data is not correct, you go
%back to the beginning to try to redo it.
answersavedata = 'No';
while matches(answersavedata, 'No')
    % Get best frame
    f = figure;
    t = tiledlayout(size(app.pwv_frames, 4),1);
    for ind = 1:size(app.pwv_frames, 4)
        nexttile;
        imagesc(app.pwv_frames(:,:,2,ind));
        title(['frame ' num2str(ind)]);
    end

    f.Position = [10 50 700 800];
    optionlist = {'Frame 1', 'Frame 2', 'Frame 3', 'Frame 4',' Frame 5'};
    [app.pwv_index, ~] = listdlg('ListString', optionlist(1:size(app.pwv_frames,4)));
    close(f)

    f = figure;
    imagesc(app.pwv_frames(:,:,2,app.pwv_index));

    clear ind optionlist t pwv_v pwv_file tmp v

    %% show red lines for ECG
    % Adjusted from change_display_index function
    start_ECG = 640;
    end_ECG = 680;
    imsize=size(app.pwv_frames(:,:,2,app.pwv_index));
    startECGline = line([1 imsize(2)],[start_ECG start_ECG],'Color','red');
    endECGline = line([1 imsize(2)],[end_ECG end_ECG],'Color','red');

    xlim([400 imsize(2)]);
    ylim([400 imsize(1)]);

    Answer = inputdlg({'ECG top', 'ECG low', 'Speed range of PW recording'},...
        'Give Info', [1 45; 1 45; 1 45], {num2str(end_ECG), num2str(start_ECG), '560'});
    close(f)
    f = figure;
    imagesc(app.pwv_frames(:,:,2,app.pwv_index));
    startECGline = line([1 imsize(2)],[str2double(Answer{2}) str2double(Answer{2})],'Color','red');
    endECGline = line([1 imsize(2)],[str2double(Answer{1}) str2double(Answer{1})],'Color','red');
    app.SpeedRangeEditField.Value = Answer{3};
    app.ECGLowEditFieldPWV.Value = str2double(Answer{2});
    app.ECGTopEditFieldPWV.Value = str2double(Answer{1});
    clear startECGline endECGline start_ECG end_ECG Answer imsize

    %% Get right part of image in order to do the PW segmentation (zoom)
    % Adjusted from SegmentPWVCurveButtonPushed
    msgbox('Zoom in on the PW curves. Press any key in the command Window when done.')
    pause
    ax = gca;
    xlimits = round(ax.XLim);
    ylimits = round(ax.YLim);
    app.Marleen.xlimits = xlimits;
    app.Marleen.ylimits = ylimits;
    cropped=squeeze(app.pwv_frames(ylimits(1):ylimits(2),xlimits(1):xlimits(2),2,app.pwv_index));
    close(f);

    % Get segmentation like in the app
    int_factor = 10;
    app.pwv_interp_frames{app.pwv_index}=imgaussfilt(imresize(cropped(:,:),int_factor),int_factor);
    vidFramecbin = bwmorph(app.pwv_interp_frames{app.pwv_index}>(mean(app.pwv_interp_frames{app.pwv_index})) & app.pwv_interp_frames{app.pwv_index}<max(app.pwv_interp_frames{app.pwv_index}),'close');

    %% Make adjustments (not in the app)
    f1 = figure;
    tiledlayout(3,1)
    nexttile
    imagesc(cropped)
    nexttile
    imagesc(vidFramecbin)
    title('Segmentation like in the app')
    % Change the imagesc to be less noisy
    vidFramecbintemp = bwareafilt(vidFramecbin, [1000 size(vidFramecbin,1)*size(vidFramecbin,2)]);
    nexttile
    imagesc(vidFramecbintemp)
    title('Segmentation with bwareafilt (1000+ pixels)')

    manualloop = 0;
    f2 = figure;
    imagesc(vidFramecbintemp)
    while manualloop == 0
        answer = questdlg('Good enough?', 'PW segmentation', 'Yes', 'No', 'No');
        close(f2)

        switch answer
            case 'Yes'
                manualloop = 1; % to exit out of while loop
                vidFramecbin = vidFramecbintemp;
                clear vidFramecbintemp
                close(f1)

            case 'No' %delete noise manually
                f2 = figure;
                imagesc(vidFramecbintemp)
                roi = drawrectangle; % get xmin, ymin, width, height
                vidFramecbintemp(roi.Position(2):roi.Position(2)+roi.Position(4), ... %y axis
                    roi.Position(1):roi.Position(1)+roi.Position(3)) = 0; % x axis
                imagesc(vidFramecbintemp)
                title('Adjusted image')
        end
    end
    clear answer ax cropped f f1 f2 int_factor manualloop roi t tmp xlimits ylimits

    %% Calculate velocity
    centroidVelocity=zeros(size(vidFramecbin,2),1);
    for j=1:size(vidFramecbin,2)
        if isempty(find(vidFramecbin(:,j),1,'last'))
            centroidVelocity(j)=1;
        else
            % This corrects for the folding effect of bad PWV
            centroidVelocity(j)=find(vidFramecbin(:,j),1,'last')-find(vidFramecbin(:,j),1,'first');
        end
    end

    %note: these numbers do not exactly fit the pw image. if you want to be
    %more precise, adjust these values. 418 - 610. It's basically the same
    %though.
    pixelSize=str2double(app.SpeedRangeEditField.Value)/((594-400)*10);  %********400 and 594 are the pixel numbers that start and end the pwv image.
    centroidVelocity=centroidVelocity*pixelSize;
    
    f = figure;
    plot(centroidVelocity)
    % Save into structure
    app.current_results.pwv_mode(app.pwv_index).velocity=centroidVelocity;
    clear j pixelSize vidFramecbin centroidVelocity

    %% Estimate ECG peaks
    % Adjusted from EstimateECGPeaksButtonPWVPushed

    % Location of ECG bar in image
    start_ECG=app.ECGLowEditFieldPWV.Value;
    end_ECG=app.ECGTopEditFieldPWV.Value;

    % xlimits = round(get(app.UIAxesPWVMode,'xlim'));
    xlimits = round(app.Marleen.xlimits);
    frame=squeeze(app.pwv_frames(:,:,2,app.pwv_index));
    % Crop image to have a region only containing ECG
    framesECG=frame(start_ECG:end_ECG,xlimits(1):xlimits(2));
    imgI = imresize(framesECG,10);
    value=zeros(size(imgI,2),1);

    for i=1:size(imgI,2)
        [pks,locs] = findpeaks(flipud(imgI(:,i)),'NPeaks',1,'MinPeakHeight',40);  %before MinPeakDistance 60 max(ECG)*0.55
        value(i)=locs(1);
    end
    value=value-min(value);
    [pks,locs] = findpeaks(value,'MinPeakHeight',max(value)*0.9);
    xline(locs);
    % hold(app.UIAxesPWVModeTmp,'on')
    % xline(app.UIAxesPWVModeTmp,locs);
    % hold(app.UIAxesPWVModeTmp,'off')
    % hold(app.UIAxesPWVModeTmp2,'on')
    % xline(app.UIAxesPWVModeTmp2,locs);
    % hold(app.UIAxesPWVModeTmp2,'off')

    app.current_results.pwv_mode(app.pwv_index).ecg_peaks=locs;
    answersavedata = questdlg('Save?', 'PW data', 'Yes', 'No', 'Yes');

    switch answersavedata
        case 'Yes'
            % Add frame to data button pushed:
            app.current_results.pwv_mode(app.pwv_index).use_data = 1;
        case 'No'
            msgbox('Redoing function, data not saved.')
    end
    
    close(f)
    clear start_ECG end_ECG framesECG frame i imgI locs pks value xlimits

end %end of while loop to see if you want to save the data or redo it.

%% Compute mean curve
% Average velocity based on consecutive cardiac cycles for each case considered
min_dim=0;
firsttime=0;
for i=1:length(app.current_results.pwv_mode)
    if(app.current_results.pwv_mode(i).use_data == 1)
        cycleLen=app.current_results.pwv_mode(i).ecg_peaks(2:end)-app.current_results.pwv_mode(i).ecg_peaks(1:end-1);
        cycleLength=min(cycleLen);
        mean_velocity=zeros(cycleLength,1);
        for j=1:length(cycleLen)
            mean_velocity=mean_velocity+app.current_results.pwv_mode(i).velocity(app.current_results.pwv_mode(i).ecg_peaks(j):app.current_results.pwv_mode(i).ecg_peaks(j)+cycleLength-1);
        end
        mean_velocity=mean_velocity/length(cycleLen);
        %diameterPixelTime = 0.0011/10;
        % diameterPixelTime = (60/min(ECGNum))/cycleLength;
        app.current_results.pwv_mode(i).mean_velocity=mean_velocity;
        if (firsttime==0)
            min_dim=length(app.current_results.pwv_mode(i).mean_velocity);
            firsttime=1;
        else
            if(length(app.current_results.pwv_mode(i).mean_velocity)<min_dim)
                min_dim=length(app.current_results.pwv_mode(i).mean_velocity);
            end
        end
    end
end

%MB: Not sure what code below changes but im assuming it's necessary...
navg=0;
mean_velocity=zeros(min_dim,1);
for i=1:length(app.current_results.pwv_mode)
    if(app.current_results.pwv_mode(i).use_data)
        mean_velocity=mean_velocity+app.current_results.pwv_mode(i).mean_velocity(1:min_dim);
        navg=navg+1;
    end
end
final_result.mean_velocity = abs(mean_velocity/navg);

% final_result=app.directoryPairs{app.currentIndex}{4};
% final_result.mean_velocity=mean_velocity;
% app.directoryPairs{app.currentIndex}{4}=final_result;
% hold(app.UIAxesPWVModeTmp2,'off')
% plot(app.UIAxesPWVModeTmp2,mean_velocity);

%% Calculate Resistence index (RI) and pulsatility index (PI)
% Adjusted from "GenerateCSVFileButton" code.
PI = (max(final_result.mean_velocity)-min(final_result.mean_velocity))/mean(final_result.mean_velocity);
RI = (max(final_result.mean_velocity)-min(final_result.mean_velocity))/max(final_result.mean_velocity);
mean_velocity = mean(final_result.mean_velocity);

disp(['mean_velocity = ' num2str(mean_velocity) ';']);
disp(['PI = ' num2str(PI) ';']);
disp(['RI = ' num2str(RI) ';']);

end


% mmode_file = "C:\Users\marle\OneDrive\Documenten\PhD\P2 - GCaMP - calcium rigidification project 1-8-22\Ultrasound Data\M32\xx - 2023-03-31-15-30-26_M32 M mode 3 left-2023-03-31-11-11-06_1.avi";
function [diameter_change, pct_diameter_change] = US_analysis_by_hand_MMode(mmode_file)
%% set up
% Adjusted from displayImages function
v=VideoReader(mmode_file);
% 3 frames
v.CurrentTime=1-1/30;
tmp=readFrame(v);
app.mmode_frames = []; %to make sure the pwv_frames is in 4D double and not uint
app.mmode_frames(1:size(tmp,1),:,:,1)=tmp;
% app.MmodeFrameIndexDropDown.Items = [app.MmodeFrameIndexDropDown.Items num2str(1)];
app.current_results.mmode(1).use_data=0;

if v.Duration >1 && v.Duration <2
    v.CurrentTime=1+1/30;
    tmp=readFrame(v);
    app.mmode_frames(:,:,:,2)=tmp;
    % app.MmodeFrameIndexDropDown.Items = [app.MmodeFrameIndexDropDown.Items num2str(2)];
    app.current_results.mmode(2).use_data=0;
end

if v.Duration>=2
    % if v.Duration>=1+20/30
    v.CurrentTime=2;
    tmp=readFrame(v);
    app.mmode_frames(1:size(tmp,1),:,:,2)=tmp;
    % app.MmodeFrameIndexDropDown.Items = [app.MmodeFrameIndexDropDown.Items num2str(2)];
    app.current_results.mmode(2).use_data=0;
end

if v.Duration>=3+1/30
    v.CurrentTime=3+1/30;
    tmp=readFrame(v);
    app.mmode_frames(1:size(tmp,1),:,:,3)=tmp;
    % app.MmodeFrameIndexDropDown.Items = [app.MmodeFrameIndexDropDown.Items num2str(3)];
    app.current_results.mmode(3).use_data=0;
end

%% start
%this while loop is to make sure that if the data is not correct, you go
%back to the beginning to try to redo it.
answersavedata = 'No';
while matches(answersavedata, 'No')
% Get best frame
f = figure;
t = tiledlayout(size(app.mmode_frames, 4),1);
for ind = 1:size(app.mmode_frames, 4)
    nexttile;
    imagesc(app.mmode_frames(:,:,2,ind));
    title(['frame ' num2str(ind)]);
end

f.Position = [10 50 700 800];
optionlist = {'Frame 1', 'Frame 2', 'Frame 3', 'Frame 4',' Frame 5'};
[app.mmode_index, ~] = listdlg('ListString', optionlist(1:size(app.mmode_frames,4)));
close(f)

f = figure;
imagesc(app.mmode_frames(:,:,2,app.mmode_index));

clear ind optionlist t mmode_file tmp v

%% show red lines for ECG
% Adjusted from change_display_index function
start_ECG = 640;
end_ECG = 680;
imsize=size(app.mmode_frames(:,:,2,app.mmode_index));
startECGline = line([1 imsize(2)],[start_ECG start_ECG],'Color','red');
endECGline = line([1 imsize(2)],[end_ECG end_ECG],'Color','red');

xlim([400 imsize(2)]);
ylim([400 imsize(1)]);

Answer = inputdlg({'ECG top', 'ECG low', 'View size'},...
    'Give Info', [1 45; 1 45; 1 45], {num2str(end_ECG), num2str(start_ECG), '1.8'});
close(f)
f = figure;
imagesc(app.mmode_frames(:,:,2,app.mmode_index));
startECGline = line([1 imsize(2)],[str2double(Answer{2}) str2double(Answer{2})],'Color','red');
endECGline = line([1 imsize(2)],[str2double(Answer{1}) str2double(Answer{1})],'Color','red');
app.ViewsizeEditField.Value = str2double(Answer{3});
app.ECGLowEditField.Value = str2double(Answer{2});
app.ECGTopEditField.Value = str2double(Answer{1});
clear startECGline endECGline start_ECG end_ECG Answer imsize

%% Get right part of image in order to do the MMode segmentation (zoom)
% Adjusted from SegmentDiameterButtonPushed
msgbox('Zoom in on the MMode curves. Press any key in the command Window when done.')
pause
ax = gca;
xlimits = round(ax.XLim);
ylimits = round(ax.YLim);
app.Marleen.xlimits = xlimits;
app.Marleen.ylimits = ylimits;
cropped=squeeze(app.mmode_frames(ylimits(1):ylimits(2),xlimits(1):xlimits(2),2,app.mmode_index));
close(f);

%make image bigger (was already in code)
int_factor = 10;
cropped=imgaussfilt(imresize(cropped(:,:),int_factor),int_factor);

cutofftop = 10; % MIGHT HAVE TO ADJUST PER ACQUISITION!!
divisionline = size(cropped,1);
cutoffbottom = 0;
answercurves = 'No, change cutoff';

while matches(answercurves, 'No, change cutoff')
    % croppedtemp = cropped<cutoff; % get low values (inside vessel = dark), set them as 1, rest as 0
    croppedtemptop = cropped(1:divisionline, :)<cutofftop;
    croppedtempbottom = cropped(divisionline+1:size(cropped,1), :)<cutoffbottom;
    croppedtemp = [croppedtemptop; croppedtempbottom];
    
    croppedtemp = bwareafilt(croppedtemp, 1); %keep only largest area (vessel)

    %check if you like the segmentation, make tiledlayout
    f1 = figure;
    tiledlayout(3,1)
    nexttile % original
    imagesc(cropped)
    nexttile % inside of vessel
    imagesc(croppedtemp)
    nexttile % inside of vessel over original
    imshowpair(croppedtemp, cropped)
    f1.Position = [10 50 500 500];

    % outline
    croppededges = edge(croppedtemp, 'canny'); %get only edges
    f2 = figure;
    imagesc(croppededges)
    croppededges = double(croppededges);
    f2.Position = [600 50 500 500];

    answercurves = questdlg('Good curves?', 'MMode data', 'Yes', 'No, retake section', 'No, change cutoff', 'Yes');

    close(f1)
    close(f2)
    clear ax f f1 f2 int_factor xlimits ylimits croppedtemptop croppedtempbottom croppedtemp
    switch answercurves
        case 'Yes'
            % Continue with the rest of the pipeline
        case 'No, retake section'
            % want to continue to the next while iteration of the big while
            % loop
        case 'No, change cutoff'
            f = figure;
            imagesc(cropped);
            pause
            % cutoff = inputdlg('New Cutoff');
            answercutoff = inputdlg({'Cutoff upper part', 'Division line', 'Cutoff lower part'});
            cutofftop = str2double(answercutoff{1});
            divisionline = str2double(answercutoff{2});
            cutoffbottom = str2double(answercutoff{3});
            close(f)
    end
end

if matches(answercurves, 'No, retake section')
    msgbox('Redoing function, data not saved.')
    clear croppededges
    continue %get to the next iteration of the while loop
end

%% Get diameter
diameter=zeros(size(croppededges,1),1);
for j=1:size(croppededges,2)
    [~,locs] = findpeaks(squeeze(croppededges(:,j)),'MinPeakDistance',size(croppededges,1)/3,'SortStr','descend');
    if size(locs,1) <= 1 % if theres no peaks and...
        if j~=1 && j~= size(croppededges,2)  %it's not first or last pixel (they have no peaks due to the edge detection)
            msgbox(["Error in segmenting diameter"; "Please make sure the top and bottom of the vessel are in the frame."])
            % disp(num2str(j))
            diameter(j) = NaN;
            continue
        else
            diameter(j) = NaN;
            continue
        end
    else
        diameter(j)=abs(locs(2)-locs(1));
    end
end

pixelSize=app.ViewsizeEditField.Value/((583-388)*10);  %********388 and 583 are the pixel numbers that start and end the mmode image. Comment MB: is actually 417-611 but comes down to the same so kept like this.
diameter=diameter*pixelSize;

% Save into structure
app.current_results.mmode(app.mmode_index).diameter=diameter;
f = figure;
plot(diameter)
clear croppededges diameter j locs pixelSize

%% Estimate ECG Peaks
% Location of ECG bar in image
start_ECG=app.ECGLowEditField.Value;
end_ECG=app.ECGTopEditField.Value;
xlimits = round(app.Marleen.xlimits);

frame=squeeze(app.mmode_frames(:,:,2,app.mmode_index));
% Crop image to have a region only containing ECG
framesECG=frame(start_ECG:end_ECG,xlimits(1):xlimits(2));
imgI = imresize(framesECG,10);
% hold(app.UIAxesMModeTmp,'on')
value=zeros(size(imgI,2),1);
for i=1:size(imgI,2)
    [pks,locs] = findpeaks(flipud(imgI(:,i)),'NPeaks',1,'MinPeakHeight',40);  %before MinPeakDistance 60 max(ECG)*0.55
    value(i)=locs(1);
end
value=value-min(value);
[pks,locs] = findpeaks(value,'MinPeakHeight',max(value)*0.9);
xline(locs);

app.current_results.mmode(app.mmode_index).ecg_peaks=locs;
answersavedata = questdlg('Save?', 'MMode data', 'Yes', 'No', 'Yes');

switch answersavedata
    case 'Yes'
        % Add frame to data button pushed:
        app.current_results.mmode(app.mmode_index).use_data = 1;
    case 'No'
        msgbox('Redoing function, data not saved.')
        %will go to the next while loop iteration in a couple of lines
end

close(f)
clear start_ECG end_ECG framesECG frame i imgI locs pks value xlimits
end % of while loop

%% Compute mean curve
% Average diameter based on consecutive cardiac cycles for each case considered
min_dim=0;
firsttime=0;
for i=1:length(app.current_results.mmode)
    if(app.current_results.mmode(i).use_data == 1)
        cycleLen=app.current_results.mmode(i).ecg_peaks(2:end)-app.current_results.mmode(i).ecg_peaks(1:end-1);
        cycleLength=min(cycleLen);
        mean_diameter=zeros(cycleLength,1);
        for j=1:length(cycleLen)
            mean_diameter=mean_diameter+app.current_results.mmode(i).diameter(app.current_results.mmode(i).ecg_peaks(j):app.current_results.mmode(i).ecg_peaks(j)+cycleLength-1);
        end
        mean_diameter=mean_diameter/length(cycleLen);
        %diameterPixelTime = 0.0011/10;
        % diameterPixelTime = (60/min(ECGNum))/cycleLength;
        app.current_results.mmode(i).mean_diameter=mean_diameter;
        if (firsttime==0)
            min_dim=length(app.current_results.mmode(i).mean_diameter);
            firsttime=1;
        else
            if(length(app.current_results.mmode(i).mean_diameter)<min_dim)
                min_dim=length(app.current_results.mmode(i).mean_diameter);
            end
        end
    end
end

%MB: Not sure what code below changes but im assuming it's necessary...
navg=0;
mean_diameter=zeros(min_dim,1);
for i=1:length(app.current_results.mmode)
    if(app.current_results.mmode(i).use_data)
        mean_diameter=mean_diameter+app.current_results.mmode(i).mean_diameter(1:min_dim);
        navg=navg+1;
    end
end
final_result.mean_diameter = mean_diameter/navg;

% final_result=app.directoryPairs{app.currentIndex}{4};
% final_result.mean_diameter=mean_diameter;
% app.directoryPairs{app.currentIndex}{4}=final_result;
% hold(app.UIAxesMModeTmp,'off')
% plot(app.UIAxesMModeTmp,mean_diameter);

%% Calculate PCT diameter change
% Adjusted from "GenerateCSVFileButton" code.
diameter_change = (max(final_result.mean_diameter)-min(final_result.mean_diameter));
pct_diameter_change = (max(final_result.mean_diameter)-min(final_result.mean_diameter))/mean(final_result.mean_diameter);

disp(['diameter_change = ' num2str(diameter_change) ';']);
disp(['pct_diameter_change = ' num2str(pct_diameter_change) ';']);

end
