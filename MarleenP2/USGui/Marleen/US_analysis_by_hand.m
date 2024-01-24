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