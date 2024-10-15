%% Example Data
% This will be the pipeline for a single imaging session. To see the
% pipeline for all imaging sessions, see
% "Preprocessing_Pipeline_CaCl_GCaMP", then "Umit_Pipeline_CaCl_GCaMP", and
% lastly, for figures "Pipeline_CaCl_GCaMP" 

% Warning: the data is pretty huge. Make sure that you have enough space
% and computing power on your computer. 

% Don't forget that besides the current files from github, you will also
% need the Labeo package Umit-master for this pipeline to work:
% github.com/LabeoTech/Umit

% Provided example data is of mouse M32, A1, R2 (a sham mouse), and M34, A1, 
% R2 (a CaCl treated mouse). 
% Change path for rawdatafolder to match where you saved the data. Same for 
% SaveFolder, but make sure you keep the structure of M../A.-R./CtxImg and 
% that your SaveFolder contains 'GCaMP' somewhere.

% M32 = Sham
DataFolder = '/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2/'; 
SaveFolder = '/home/mbakker/P2_scripts/MarleenP2/ExampleData/GCaMP/M32/A1-R2/CtxImg/';

% %M34 = CaCl
% DataFolder = '/home/mbakker/P2_scripts/MarleenP2/Example/ExampleData/M34-A1-R2/'; 
% SaveFolder = '/home/mbakker/P2_scripts/MarleenP2/Example/ExampleData/GCaMP/M34/A1-R2/CtxImg/';

mkdir(SaveFolder);

% Note: in the pipeline, it will check if you already did a function. Here,
% it will not check for that, so be sure that you don't do functions twice!



%% Preprocessing_Pipeline_CaCl_GCaMP:
%% Sort colours
ImagesClassification(DataFolder, SaveFolder, 1, 1, 0);
% This will give you a files called 'fluo_475.dat', 'fluo_567.dat', 
% 'green.dat', 'red.dat', all corresponding matfiles, and 'AcqInfos.mat'. 
% The only reason that fluo_475.dat is there is because of the way the
% system saves data. This data is nonsense though, as our channel is
% fluo_567. You can delete the unnecessary data but ONLY AFTER DOING THE
% DICHROIC CORRECTION. Otherwise, code will bug. 

% You can check the data like this:
fid = fopen([SaveFolder 'fluo_567.dat']); %change if you want to see green or red
dat = fread(fid, inf, '*single' );
dat = reshape(dat, 512, 512, []); %these acquistions were 512*512 pixels
imagesc(dat(:,:,23)); % last number is frame

% for a little movie:
figure;
for ind = 1:200 
    imagesc(dat(:,:,ind), [0 40000])
    title(['frame ' num2str(ind)])
    pause(0.05)
end

clear dat

%% Mark movement based on treadmill
MarkMovedFrames(DataFolder, SaveFolder);
% This will give you the variable "MovMask.mat". If you open it, you can
% see when the mouse had detected movement on the threadmill:
load([SaveFolder 'MovMask.mat'])
figure
plot(MovMask)
% 0 means movement, 1 means no movement.

clear MovMask

%% Correct for shift between colours because of dichroic mirror
% This will do the coregistration between different colour channels, to
% make sure the brain is at the same place in the image. It will give you
% two images, and then it will pause and say in the command window "show
% video". The code will pause until you PRESS ENTER IN THE COMMAND WINDOW,
% after which a video jumping from one image to another will pop up. This
% will take a couple of seconds, and at the end a prompt will come up
% asking of the coregistration was correct. If not, you can coregister the
% images yourself. 
% The red.dat file will be overwritten with the coregistered data (as will
% the fluo_475 but we don't care about that one). 
CorrectDichroic(SaveFolder)

% see sort colours
delete([SaveFolder 'fluo_475.dat'])
delete([SaveFolder 'fluo_475.mat'])

%% flip left/right 
% Our imaging system mirrors the images, so we have to flip them to have
% the left hemisphere on the left and the right on the right. This
% overwrites the .dat files. MAKE SURE YOU ONLY DO THIS ONCE. 

for datfile = {'green', 'red', 'fluo_567'}
    fid = fopen([SaveFolder char(datfile) '.dat']);
    dat = fread(fid, inf, '*single');
    dat = reshape(dat, 512, 512, []);
    dat = fliplr(dat);

    disp(['Overwriting ' char(datfile) '.dat file'])
    fclose(fid);
    fid = fopen([SaveFolder char(datfile) '.dat'], 'w');
    fwrite(fid,dat,'single');
    fclose(fid);
end

% This is nice data for Figure 1 F. You flipped fluo last so you should
% have that one loaded. We will plot the ROI over that later.
Anatomical_image = mean(dat,3);
figure
imagesc(Anatomical_image)
save([SaveFolder 'Anatomical.mat'], 'Anatomical_image')

clear dat

%% Make ROI and mask
fid = fopen([SaveFolder 'fluo_567.dat']);
fluo_im = fread(fid, 512*512, '*single' );
fluo_im = reshape(fluo_im, 512, 512);
fclose(fid);

fid = fopen([SaveFolder 'green.dat']);
green_im = fread(fid, 512*512, '*single' );
green_im = reshape(green_im, 512, 512);
imagesc(green_im)
fclose(fid);

%increase contrast
green_im = green_im./mean(green_im(:));
green_im = (green_im-min(green_im(:)))./(max(green_im(:))-min(green_im(:)));
green_im(green_im < 0) = 0;
green_im(green_im > 1) = 1;
green_im = adapthisteq(green_im);
figure
imagesc(green_im)

% Follow instructions of the ROImanager (it will give a pop-up window with
% steps of what to do. Note: Save in the folder of the mouse, NOT in the
% CtxImg or A1-R2. Press enter in the command window if you're done. 
ROImanager(fluo_im)
f = msgbox(["In ROImanager, do the following steps:"; ...
    "-  Image – Set origin – New – drag to bregma, right mouse-click to set"; ...
    "-  Image – Set origin – Align image to origin – drag to lambda, right mouse-click to set and correct for tilt in frame"; ...
    "-  Image – Set pixel size – 50 pixels per mm"; ...
    "-  Image – Mask – Draw new – Draw the mask by clicking points, you can still drag the points after setting them, double click inside the mask to confirm"; ...
    "-  Image – Image reference file… - Export – save as ImagingReferenceFrame.mat in the mouse folder"; ...
    "-  Create – Mouse Allen Brain Atlas – Areas – Select areas – close window to select them all – drag atlas to what seems to be fitting – double click inside atlas to confirm"; ...
    "-  File – Save as… - ROImasks_data.mat in mouse folder"]);
input('Make reference and ROI')

% The allen atlas has many small regions, but since we don't do any
% verification based on stimulus data, we can never be that sure about the
% placement of these regions. Therefore, we pool regions into larger ROI. 
ClusterRois(SaveFolder, 1)
% If you rotated the brain a bit to be straight, this will make sure your
% ROI fit on your data:
Correct_for_rotation_ROI(SaveFolder)

%% Figure 1F
% In this case, you saved the anatomical image earlier. In the pipeline, we
% don't do that. However, you can always make an anatomical image by
% opening the fluo_567.dat data, since we rename the working data to
% hemoCorr_fluo.dat after this step. 
seps = strfind(SaveFolder, filesep);
load([SaveFolder(1:seps(end-2)) 'BigROI.mat'], 'AtlasMask')
load([SaveFolder 'Anatomical.mat'], 'Anatomical_image') 
AtlasMask(AtlasMask~=0) = 1;
axis('image'); axis('equal');

saveas(gcf, [SaveFolder 'Anatomical_image_with_atlas.svg'], 'svg');



%% Umit_Pipeline_CaCl_GCaMP
% Is called Umit because it's based on LabeoTech's Umit-toolbox: https://github.com/LabeoTech/Umit

%% Hemodynamic correction
% Be careful to only do this once! You'll get the file hemoCorr_fluo.dat
% out of this (with corresponding matfile)
fid = fopen([SaveFolder 'fluo_567.dat']);
dat = fread(fid, inf, '*single');
metadata = load([SaveFolder 'fluo_567.mat']);
hemoCorr_fluo = HemoCorrection(SaveFolder, dat, metadata, {'Red','Green'});
fclose(fid);

hemoCorr_fluo(isnan(hemoCorr_fluo)) = 0;

fid = fopen([SaveFolder 'hemoCorr_fluo.dat'],'w');
fwrite(fid, hemoCorr_fluo, '*single');
fclose(fid);

Make_HemoCorr_Matfile(SaveFolder);

clear dat hemoCorr_fluo

%% Save without filtering
% This is to verify that filtering doesn't do weird things. The summation
% of all mice for this (triggered on gcamp) is in appendix 3B. 
fid = fopen([SaveFolder 'hemoCorr_fluo.dat']);
dat = fread(fid, inf, '*single');
dat(dat == 0) = NaN;
dat = reshape(dat, 512, 512, []);
fclose(fid);
% normalize so you can compare to filtered later
datnofilt = dat./mean(dat,3,'omitnan');

fid = fopen([SaveFolder 'fluo_nofilt.dat'],'w');
fwrite(fid,datnofilt,'*single');
fclose(fid);

%% Normalisation & filtering of the data
hemoCorr_fluo = NormalisationFiltering(SaveFolder, 'hemoCorr_fluo', 0.3, 3, 1, 0);

fid = fopen([SaveFolder 'hemoCorr_fluo.dat'],'w');
fwrite(fid,hemoCorr_fluo,'*single');
fclose(fid);

% Visualise a bit (this will differ from appendix 3B because it's not based
% on triggered data).
% Note: this is the same pixel and timing as figure 1D, but you can change
% to explore a bit, just make sure you are inside the brain.
coords = [300 300];
timing = [800 1100];
figure
plot(squeeze(hemoCorr_fluo(coords(1), coords(2), timing(1):timing(2))));
hold on
plot(squeeze(datnofilt(coords(1), coords(2), timing(1):timing(2))));
legend({'filtered', 'unfiltered'})

clear dat datnofilt hemoCorr_fluo

%% Calculate HbO and HbR
% This will save as .dat files again
[~, ~] = HemoCompute(SaveFolder, SaveFolder, 'gcamp', {'red', 'green'}, 1);

%% Figure 1D - if you plot for M32
% With this data you can replicate Figure 1D:
% (For cleaner way of plotting, see ExampleData_GCaMP_HbO_HbR, this is just 
% to visualize quickly)
coords = [300 300];
timing = [800 1100];

% (Open hemoCorr_fluo if you dont have it still open from function above)
fid = fopen([SaveFolder 'hemoCorr_fluo.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);
fclose(fid);

figure
yyaxis left
plot(squeeze(dat(coords(1), coords(2), timing(1):timing(2))));

yyaxis right
fid = fopen([SaveFolder 'HbO.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);
fclose(fid);
plot(squeeze(dat(coords(1), coords(2), timing(1):timing(2))));
hold on
fid = fopen([SaveFolder 'HbR.dat']);
dat = fread(fid, inf, '*single');
dat = reshape(dat, 512, 512, []);
fclose(fid);
plot(squeeze(dat(coords(1), coords(2), timing(1):timing(2))));
clear dat
legend('gcamp', 'hbo', 'hbr')

%% Detect outliers
% Detect extreme frames, label them as outliers. Can be because of movement
% for example. Do the same for pixels. With pixels, it can be that the back
% of the mouse got into frame often at the bottom of the brain, and thus
% those pixels are better excluded. 
MakeOutlierMask(SaveFolder, {'hemoCorr_fluo', 'HbO', 'HbR'}, 0);

% Visualize:
load([SaveFolder 'OutlierMask.mat'])
fid = fopen([SaveFolder 'hemoCorr_fluo.dat']);
dat = fread(fid, 512*512, '*single');
dat = reshape(dat, 512, 512, []);
fclose(fid);

figure
tiledlayout('flow');
nexttile
imagesc(dat, [0.95 1.05])
title('first GCaMP frame')
nexttile
imagesc(OutlierPixels.HbO)
title('HbO')
nexttile
imagesc(OutlierPixels.HbR)
title('HbR')
nexttile
imagesc(OutlierPixels.hemoCorr_fluo)
title('GCaMP')
nexttile
imagesc(dat.*OutlierPixels.HbO.*OutlierPixels.HbR.*OutlierPixels.hemoCorr_fluo, [0.95 1.05])
title('Outliers Removed')

%% Get timecourses
% Timecourses are used for the correlation matrices. It's based on a seed
% in the region of interest, and can later be correlated with the other
% seeds to get the functional connectivity. 

% The manual input is 1. It will give you the suggested middle of the
% regions of interest of the allen atlas. However, sometimes this is over a
% bubble or damage in the window, which would throw off the data. You can
% therefore alter the x and y coordinates of certain seeds. It will ask to
% confirm at the end. If the code seems stuck, hit a key in the command
% window. It will save the seeds as Seeds.mat, and the timecourses as
% timecourses_hemoCorr_fluo_centroids.mat or similar.
GetTimecourses(SaveFolder, 'hemoCorr_fluo', 1)
GetTimecourses(SaveFolder, 'HbO', 1)
GetTimecourses(SaveFolder, 'HbR', 1)

% With GSR:
GetTimecourses(SaveFolder, 'hemoCorr_fluo', 1, 1)
GetTimecourses(SaveFolder, 'HbO', 1, 1)
GetTimecourses(SaveFolder, 'HbR', 1, 1)



%% Figure 2, 3A, and 4A
% Get average curves. Normally this function is inside other functions (for
% example plotNVCinRS). This should give you NVC_ROI and NVC_ROI_random. In
% NVC_ROI are the curves around detected activations
GetAverageCurves(SaveFolder)

% To replicate Figure 2 (or 4A) with a single mouse, put WholeBrain as region. 
% For Figure 3, change the regions to what you prefer:
load([SaveFolder 'NVC_ROI.mat'])
region = 'VisualROI_R'; % can be any other region
figure
title(region, 'Interpreter', 'none')
yyaxis left
plot(fluocurves.(region), 'color', 'green')
yyaxis right
plot(hbocurves.(region), 'color', 'red')
hold on
plot(hbrcurves.(region), 'color', 'blue', 'LineStyle', '-')

% If you want, you can compare it to the NVC_ROI_random (not depicted in
% paper, only discussed in text). 
load([SaveFolder 'NVC_ROI_random.mat'])
figure
title([region ' random'], 'Interpreter', 'none')
yyaxis left
plot(fluocurves_random.(region), 'color', 'green')
yyaxis right
plot(hbocurves_random.(region), 'color', 'red')
hold on
plot(hbrcurves_random.(region), 'color', 'blue', 'LineStyle', '-')

%% Figure 4B
% Get spontaneous activations based on HbO. If you want the "random"
% detected activations, you can do:
% GetAverageCurves(SaveFolder, 'hemoCorr_fluo', 'NVC_ROI_HbO', 'HbO', 'HbO', 'HbR', 1)
% However, this will give similar results as the random GCaMP activations,
% with the only difference being the number of activations (to match their
% respective datatypes).
GetAverageCurves(SaveFolder, 'hemoCorr_fluo', 'NVC_ROI_HbO', 'HbO')

load([SaveFolder 'NVC_ROI_HbO.mat'])
region = 'WholeBrain'; % can be any other region
figure
title([region ' - HbO-detected'], 'Interpreter', 'none')
yyaxis left
plot(fluocurves.(region), 'color', 'green')
yyaxis right
plot(hbocurves.(region), 'color', 'red')
hold on
plot(hbrcurves.(region), 'color', 'blue', 'LineStyle', '-')

%% Appendix 3B
ExampleData_NoFilt(SaveFolder, SaveFolder)
% it will save a .svg figure in a folder called SingleSubjectExample

%% Figure 6A
% Correlation matrices
% These figures will differ from the ones in the paper because you are now
% plotting for only one acquisition, whereas the one in the paper is the
% average of all mice (in a certain group).

load('/home/mbakker/P2_scripts/MarleenP2/NL.mat'); % CHANGE to fit your own path

datatypes = {'timecourses_hemoCorr_fluo_centroids.mat', ...
    'timecourses_HbO_centroids.mat', 'timecourses_HbR_centroids.mat',...
    'timecourses_hemoCorr_fluo_centroids_GSR.mat', ...
    'timecourses_HbO_centroids_GSR.mat', 'timecourses_HbR_centroids_GSR.mat'};

% datatypes = {'timecourses_hemoCorr_fluo_centroids.mat', ...
%     'timecourses_hemoCorr_fluo_centroids_GSR.mat'};

% datatypes = {'timecourses_hemoCorr_fluo_centroids.mat', 'timecourses_hemoCorr_fluo_centroids_GSR.mat'};
for ind = 1:size(datatypes,2)
    figure
    load(datatypes{ind})
    Timecourses = cell2mat(AllRois(:,2));
    load([SaveFolder 'MovMask.mat'], 'MovMask');
    Timecourses = Timecourses .* MovMask; 
    load([SaveFolder 'OutlierMask.mat'], 'OutlierFrames');
    Timecourses = Timecourses .* OutlierFrames;
    Timecourses(Timecourses == 0) = nan;
    Timecourses = Timecourses(:,1:9000);
    Timecourses(3,:) = []; % remove auditory cortex (often outside window)
    Timecourses(7,:) = [];

    data = corr(Timecourses', 'rows','pairwise');

    % all this below is to make it the same order and colour as in the paper,
    % but you could also just do: imagesc(data)
    allrois = {'VR', 'SR', 'MR', 'RR', 'VL', 'SL', 'ML', 'RL'};
    ax = gca;
    [allrois, sortindex] = sort(allrois);
    data = reshape(data, size(allrois,2), size(allrois,2));
    data = data(sortindex, sortindex);
    data = tril(data);
    data(data == 0 ) = NaN;
    data(data == 1 ) = NaN;

    if contains(datatypes{ind}, 'GSR')
        imagesc(ax, data, 'AlphaData', ~isnan(data), [-1 1])
    else
        imagesc(ax, data, 'AlphaData', ~isnan(data), [0 1])
    end
    ax.Box = 'off';
    axis image
    axis equal;

    yticks(1:size(allrois,2));
    yticklabels(allrois);
    ay = get(gca,'YTickLabel');
    set(gca, 'YTickLabel', ay);
    xticks(1:size(allrois,2));
    xticklabels(allrois);
    xtickangle(90)

    colormap(NL)

    title(datatypes{ind}, 'Interpreter', 'none')
end

% That's it :)
% For group plots, calculations or statistics, see the actual pipelines
% (Pipeline_CaCl_GCaMP).

