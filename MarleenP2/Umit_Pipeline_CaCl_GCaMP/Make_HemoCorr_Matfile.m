% in a seperate function so you don have all these variables in your
% workspace.

function Make_HemoCorr_Matfile(DataFolder)

% have to save a .mat file for this particular .dat file
load([DataFolder 'fluo_567.mat']);
datFile = 'hemoCorr_fluo.dat';

if ~exist('Stim', 'var')
    Stim = 0;
end

save([DataFolder 'hemoCorr_fluo.mat'], 'Datatype', 'FirstDim', ...
    'Freq', 'Stim', 'datFile', 'datLength', 'datName', 'datSize', ...
    'dim_names', 'tExposure', '-v7.3');
end