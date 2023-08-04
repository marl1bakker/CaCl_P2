function CompareMovementRecordings(directory)

if( ~strcmp(directory(end), filesep) )
    directory = [directory filesep];
end

%% Make acquisition list
%get list of acquisitions
[Mice, AcqList] = MakeAcqList(directory);

if ~isfile('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
    %make table to store the right recordings in the end
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    A1 = cellstr(repmat('empty', size(Mice,1),1));
    A2 = cellstr(repmat('empty', size(Mice,1),1));
    A3 = cellstr(repmat('empty', size(Mice,1),1));
    Sex = cellstr(repmat('empty', size(Mice,1),1));
    Group = cellstr(repmat('empty', size(Mice,1),1));
    
    RecordingOverview = table(Mouse, Folder, A1, A2, A3, Sex, Group);
    clear Mouse A1 A2 A3 Folder Sex Group
    
    
else
    %make sure you add mice if you have more, but don't overwrite old ones
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
    A1 = repmat('empty', 1,1);
    A2 = repmat('empty', 1,1);
    A3 = repmat('empty', 1,1);
    Sex = repmat('empty', 1,1);
    Group = repmat('empty', 1,1);
    
    for ind = 1:size(Mouse, 1)
        if ~matches(Mouse(ind), 'M1') && ~matches(Mouse(ind), 'M2') %M1 and M2 are special, skip for now
            %search if already done
            if sum(matches(RecordingOverview.Mouse, Mouse(ind))) > 0
                disp([Mouse(ind) ' already in table'])
            else %if not, add mouse to end
                %             Mouse = Mouse(ind);
                %             Folder = Folder(ind);
                AddMouse = {Mouse(ind), Folder(ind), A1, A2, A3, Sex, Group};
                RecordingOverview = [RecordingOverview; AddMouse];
            end
        end
    end
end

clear AddMouse A1 A2 A3 Folder Group ind Mouse Sex Mice

%% Add group and sex for mice
load('/home/mbakker/P2_scripts/MarleenP2/MiceCodes.mat');

for ind = 1:size(RecordingOverview.Mouse,1)
    MouseIndex = find(strcmp(Mice.CodeOfMouse, RecordingOverview.Mouse(ind)));
    
    if isequal(RecordingOverview.Group(ind), {'empty'})
        RecordingOverview.Group(ind) = Mice.CaClSham(MouseIndex);
        RecordingOverview.Sex(ind) = Mice.MaleFemale(MouseIndex);
    end
end

%% Start movement checking
%these lists are only for the current directory that was given, and not the
%whole recordingoverview!
AcqListMice = extractBetween(struct2cell(AcqList(:)), '/GCaMP/', '/')';
AcqListAcq = extractBetween(struct2cell(AcqList(:)), AcqListMice', 'R')';
AcqListRec = extractAfter(struct2cell(AcqList(:)), AcqListAcq')';
AcqListAcq = extractBetween(AcqListAcq, '-', '-');

% AcqList = reshape(struct2cell(AcqList), 1, []);

Acquisitions = {'A1','A2','A3'};


for ind1 = 1:size(RecordingOverview.Mouse) %go per mouse
    if isequal(RecordingOverview.A3(ind1), {'empty'}) %assume that if A3 is done, all are done. Don't redo them
        mouse = RecordingOverview.Mouse{ind1};
        savefolder = RecordingOverview.Folder{ind1};
        disp(mouse)
        
        for ind2 = 1:size(Acquisitions,2) %go per acquisition
            acq = Acquisitions{ind2};
            disp(acq)
            %             Recordings = {};
            
            %find all recordings of that acquisition
            indexforacqlist = strcmp(mouse, AcqListMice) .* strcmp(acq, AcqListAcq);
            indexforacqlist = find(indexforacqlist);
            
            for ind3 = 1:size(indexforacqlist)
                savename = ['-' mouse '-' mouse '-' acq '-' char(AcqListRec(indexforacqlist(ind3)))];
                %                 Recordings{end+1} = AcqList(indexforacqlist(ind3)).name;
                
                if exist([savefolder filesep mouse filesep 'Movement' savename '.tif'], 'file')
                    % if we already saved it, just show that one
                    MovIm = imread([savefolder filesep mouse filesep 'Movement' savename '.tif']);
                    figure('Position', [50 200 400 400]);
                    imagesc(MovIm)
                else
                    % otherwise, plot the movement
                    %                     PlotTreadmill(AcqList(index).name) %plot R1, R2, R3 etc.
                    PlotTreadmill(AcqList(indexforacqlist(ind3)).name) % plot R1, R2, R3 etc.
                    saveas(gcf, [savefolder filesep mouse filesep 'Movement' savename], 'tiff');
                end
                
            end
            Recordings = struct2cell(AcqList(indexforacqlist));
            [indx, ~] = listdlg('PromptString', {'Which recording do you want to use'},...
                'SelectionMode','single','ListString', Recordings, 'ListSize', [600 300]);
            
            % add right recording to the overview
            eval(['RecordingOverview.' acq ' (ind1) = Recordings(indx);']);
            close all
        end
    end
end

save( '/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview')

end











% bla = dir('GCaMPMovement*');
% 
% for ind = 1:size(bla,1)
%     newfilename = extractAfter(bla(ind).name, 'GCaMP');
%     mouse = extractBetween(bla(ind).name, 'GCaMPMovement-', '-');
%     movefile(bla(ind).name, ['/media/mbakker/PJM - HDD - 2/Marleen/GCaMP/' char(mouse) filesep newfilename])
% end







%% If you want to mark frames for movement, check Pipeline_CaCl_GCaMP or do:
%
% Treadmill = (Treadmill < 1.6455) & (Treadmill > 1.64); % get a logical array for movement (0) or not (1)
%
% load([DataFolder 'green.mat']);
% TimestepsPerFrame = size(Treadmill,1)/datLength;
%
% FramesMoved = ones(datLength,1);
% ind1 = 1; %just to get it started
% for ind = 1:datLength
%     ind2 = ceil(ind * TimestepsPerFrame);
%
%     if any(Treadmill(ind1:ind2)<1) % gives 1 if there was movement
%         FramesMoved(ind) = 0; %Give a 0 for frames that have movement
%     end
%
%      ind1 = floor(ind * TimestepsPerFrame) ; %this is already updated for next step
% end
