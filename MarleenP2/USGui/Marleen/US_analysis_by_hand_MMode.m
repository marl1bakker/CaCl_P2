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