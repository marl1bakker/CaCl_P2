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