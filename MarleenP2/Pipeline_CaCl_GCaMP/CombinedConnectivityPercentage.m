function CombinedConnectivityPercentage(dataname, Acquisition, Grouping)

%% set up
if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo';
end

if ~exist('Acquisition', 'var')
    Acquisition = 'A1';
end

if ~exist('dataname', 'var')
    dataname = 'hemoCorr_fluo.dat';
elseif  length(dataname) < 4 || ( ~strcmp(dataname(end-3:end), '.dat') )
    dataname = [dataname '.dat'];
end

load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview');


%% Choose grouping
if ~exist('Grouping', 'var')
    Possibilities = {'Group', 'Sex'};
    [indx, ~] = listdlg('ListString', Possibilities);
    Grouping = Possibilities(indx);
end

[groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping);

if size(Grouping, 2)>1
    Grouping = 'Combi';
else
    Grouping = char(Grouping);
end
%% Single group Corr matrix
% Start going per group, per mouse
for index = 1:size(groups,1) % Go per group
    
    group = groups{index};
    disp(group);
    eval(['idx = RecordingOverview.' Grouping ' == ''' group ''';'])
    Mousegroup = RecordingOverview(idx,:);
    eval([group 'N = 0;']) %this is not always equal to the Mousegroup size, because maybe we don't have the Timecourse for that mouse yet. That's why the n group is seperately calculated.
    AllConnValues = [];
    AllPercentages = [];
    
    if exist(['/media/mbakker/GDrive/P2/GCaMP/ConnectivityPercentage/AllConnValues_' group '.mat'], 'file')
        load(['/media/mbakker/GDrive/P2/GCaMP/ConnectivityPercentage/AllConnValues_' group '.mat'], 'AllConnValues', 'AllPercentages');
        
        eval(['AllConnValues_' group '= AllConnValues;']);
        eval(['AllPercentages_' group '= AllPercentages;']);
        
    else
        
        for ind = 1:size(Mousegroup, 1) %go per mouse
            Mouse = Mousegroup.Mouse{ind};
            eval(['DataFolder = [Mousegroup.' Acquisition '{ind} filesep];']);
            SaveFolder = [Mousegroup.SaveDirectory{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
            
            [ConnValuesMouse, nrpixels] = SingleSubjectConnectivityValues(SaveFolder, dataname);
            PercentagesMouse = ConnValuesMouse/nrpixels;
            
            AllConnValues (ind,:) = ConnValuesMouse;
            AllPercentages(ind,:) = PercentagesMouse;
        end
        
        save(['/media/mbakker/GDrive/P2/GCaMP/ConnectivityPercentage/AllConnValues_' group], 'AllConnValues', 'AllPercentages');
        
        eval(['AllConnValues_' group '= AllConnValues;']);
        eval(['AllPercentages_' group '= AllPercentages;']);
        
    end
end

clear AllConnValues AllPercentages idx index indx Possibilities 

end