%         RecordingOverview.Combi = cellstr(repmat(' ', size(RecordingOverview,1),1));
%         for ind = 1:size(Grouping, 2)
%             eval(['currentgroup = cellstr(RecordingOverview.' Grouping{ind} ');' ]);
%             RecordingOverview.Combi = strcat(RecordingOverview.Combi, currentgroup);
%         end
%         groups = unique(RecordingOverview.Combi);
%         Grouping = 'Combi';
%         RecordingOverview.Combi = categorical(RecordingOverview.Combi);

function CompareMovementRecordings(directory, imagingtype)
if ~exist('imagingtype', 'var')
    imagingtype = 'GCaMP';
end

if( ~strcmp(directory(end), filesep) )
    directory = [directory filesep];
end

if matches(imagingtype, 'GCaMP')
    savenamerecordingoverview = '/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat';
else
    savenamerecordingoverview = ['/home/mbakker/P2_scripts/MarleenP2/RecordingOverview_' imagingtype '.mat'];
end

%% Add mice to Recordingoverview
%get list of acquisitions
[Mice, AcqList] = MakeAcqList(directory);

if ~isfile(savenamerecordingoverview)
    %make table to store the right recordings in the end
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    A1 = cellstr(repmat('empty', size(Mice,1),1));
    A2 = cellstr(repmat('empty', size(Mice,1),1));
    A3 = cellstr(repmat('empty', size(Mice,1),1));
    Group = categorical(cellstr(repmat('empty', size(Mice,1),1)));
    Sex = categorical(cellstr(repmat('empty', size(Mice,1),1)));
    Combi = categorical(cellstr(repmat('empty', size(Mice,1),1)));
    SaveDirectory = cellstr(repmat('empty', size(Mice,1),1));
    
    RecordingOverview = table(Mouse, Folder, A1, A2, A3, Group, Sex, Combi, SaveDirectory);
    clear Mouse A1 A2 A3 Folder Sex Group SaveDirectory
    
else
    %make sure you add mice if you have more, but don't overwrite old ones
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    load(savenamerecordingoverview, 'RecordingOverview')
    A1 = repmat('empty', 1,1);
    A2 = repmat('empty', 1,1);
    A3 = repmat('empty', 1,1);
    Group = categorical(cellstr(repmat('empty', 1,1)));
    Sex = categorical(cellstr(repmat('empty', 1,1)));
    Combi = categorical(cellstr(repmat('empty', 1,1)));
    SaveDirectory = repmat('empty', 1,1);
    
    for ind = 1:size(Mouse, 1)
        if matches(Mouse(ind), 'M1') || matches(Mouse(ind), 'M2')
            continue
        elseif sum(matches(RecordingOverview.Mouse, Mouse(ind))) > 0
%             disp([Mouse{ind} ' already in table'])
            continue
        else
            AddMouse = {Mouse(ind), Folder(ind), A1, A2, A3, Group, Sex, Combi, SaveDirectory};
            RecordingOverview = [RecordingOverview; AddMouse];
        end
    end
end
clear AddMouse A1 A2 A3 Folder Sex Group SaveDirectory ind Mouse Sex Mice

%% Add group and sex for mice
load('/home/mbakker/P2_scripts/MarleenP2/MiceCodes.mat', 'Mice');

for ind = 1:size(RecordingOverview.Mouse,1)
    MouseIndex = find(strcmp(Mice.CodeOfMouse, RecordingOverview.Mouse(ind)));
    
    if isequal(RecordingOverview.Group(ind), {'empty'})
        RecordingOverview.Group(ind) = Mice.CaClSham(MouseIndex);
        RecordingOverview.Sex(ind) = Mice.MaleFemale(MouseIndex);
        combiname = strcat(cellstr(RecordingOverview.Group(ind)), cellstr(RecordingOverview.Sex(ind)));
        RecordingOverview.Combi(ind) = categorical(combiname);
    end
end

%% Add Saving directory
AcqListMice = extractBetween(struct2cell(AcqList(:)), ['/' imagingtype '/'], '/')';
AcqListAcq = extractBetween(struct2cell(AcqList(:)), AcqListMice', 'R')';
AcqListRec = extractAfter(struct2cell(AcqList(:)), AcqListAcq')';
AcqListAcq = extractBetween(AcqListAcq, '-', '-');

Acquisitions = {'A1','A2','A3'};

for ind = 1:size(RecordingOverview.Mouse, 1) % go per mouse
    mouse = RecordingOverview.Mouse{ind};
    savefolder = RecordingOverview.Folder{ind};
%     disp(mouse)
    
    if matches(mouse, 'M1') || matches(mouse, 'M2') % mouse 1 and 2 are special, skip for now
        continue
    end
    
    % add save directory
    Mousenr = str2num(mouse(2:end));
    if (Mousenr >= 15 && Mousenr <= 36) && matches(RecordingOverview.SaveDirectory{ind}, 'empty')
        RecordingOverview.SaveDirectory{ind} = ['/media/mbakker/GDrive2/P2/' imagingtype];
    elseif matches(RecordingOverview.SaveDirectory{ind}, 'empty')
        RecordingOverview.SaveDirectory{ind} = ['/media/mbakker/GDrive/P2/' imagingtype];
    end
        
    for ind2 = 1:size(Acquisitions,2) %go per acquisition
        acq = Acquisitions{ind2};
        if eval(['~matches(RecordingOverview.' acq '(ind), ''empty'');'])
%             disp(['Acquisition ' mouse ' ' acq ' already done'])
            continue
        end
            
        %find all recordings of that acquisition
        indexforacqlist = strcmp(mouse, AcqListMice) .* strcmp(acq, AcqListAcq);
        indexforacqlist = find(indexforacqlist);
            
        %get plots
        for ind3 = 1:size(indexforacqlist)
            savename = ['-' mouse '-' mouse '-' acq '-' char(AcqListRec(indexforacqlist(ind3)))];
            
            if exist([savefolder filesep mouse filesep 'Movement' savename '.fig'], 'file')
                open([savefolder filesep mouse filesep 'Movement' savename '.fig'])
            else
                PlotTreadmill(AcqList(indexforacqlist(ind3)).name) 
                saveas(gcf, [savefolder filesep mouse filesep 'Movement' savename '.fig']);
            end
        end
        
        Recordings = struct2cell(AcqList(indexforacqlist));
        
        if size(Recordings, 3)>1
            [indx, ~] = listdlg('PromptString', {'Which recording do you want to use'},...
                'SelectionMode','single','ListString', Recordings, 'ListSize', [600 300]);
            eval(['RecordingOverview.' acq ' (ind) = Recordings(indx);']);
        elseif size(Recordings,3) == 1
            eval(['RecordingOverview.' acq ' (ind) = Recordings(1);']);
        end
        
        close all
    end
end

save( savenamerecordingoverview, 'RecordingOverview')

end

function [Mice, AcqList] = MakeAcqList(directory)

if ~exist('directory', 'var')
    directory = '/media/mbakker/data1/P2/GCaMP';
end

if( ~strcmp(directory(end), filesep) )
    directory = [directory filesep];
end

% make a list of folders with acquisitions and mice
NrOfMice = size(dir([directory 'M*']),1);
Mice = dir([directory 'M*']);
% Mice = Mice(3:end);
AcqList = {};


for index = 1:NrOfMice
    Mouse = Mice(index).name;
    
    DataFolder = [directory Mouse filesep];
    AcqMouse = dir([DataFolder 'M*']); % get all the folders/acquisitions of a certain mouse
    AcqMouse = AcqMouse([AcqMouse.isdir] == 1); % get only folders
    
    for ind = 1:(size(AcqMouse,1)) %
        AcqList{end+1} = [DataFolder AcqMouse(ind).name];
    end
    
end

% Mice = struct('name', Mice);
AcqList = struct('name', AcqList);

end

function PlotTreadmill(DataFolder)

UpperLimit = 1.65;
LowerLimit = 1.63;

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

aiFilesList = dir([DataFolder 'ai_*.bin']);

% Opening of the files:
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

Treadmill = AnalogIN(:,5);

x = linspace(0, 600, size(Treadmill,1));

figure()
plot(x', Treadmill)
ylim([1.62 1.66])
title(DataFolder(end-9:end-1))
%to set lines, but the range changes per acq
% line([1 size(Treadmill,1)], [UpperLimit UpperLimit], 'Color', 'red') 
% line([1 size(Treadmill,1)], [LowerLimit LowerLimit], 'Color', 'red')

NotMoved = sum((Treadmill < UpperLimit) & (Treadmill > LowerLimit)); % get a logical array for movement (0) or not (1)
Moved = size(AnalogIN,1) - NotMoved;
FractionMoved = Moved/size(AnalogIN,1);

subtitle(['Fraction movement = ' num2str(FractionMoved)])

end
