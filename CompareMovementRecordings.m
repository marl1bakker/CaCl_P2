function CompareMovementRecordings(directory)

if( ~strcmp(directory(end), filesep) )
    directory = [directory filesep];
end

%get list of acquisitions
[Mice, AcqList] = MakeAcqList(directory);

if ~isfile('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
    %make table to store the right recordings in the end
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    A1 = cellstr(repmat('empty', size(Mice,1),1));
    A2 = cellstr(repmat('empty', size(Mice,1),1));
    A3 = cellstr(repmat('empty', size(Mice,1),1));

    RecordingOverview = table(Mouse, Folder, A1, A2, A3);
    clear Mouse A1 A2 A3 Folder
    
else 
    %make sure you add mice if you have more, but don't overwrite old ones
    Mouse = struct2table(Mice).name;
    Folder = struct2table(Mice).folder;
    load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
    for ind = 1:size(Mouse, 1)
        %search if already done
        if sum(contains(RecordingOverview.Mouse, Mouse(ind))) >= 0
            disp([Mouse(ind) ' already in table'])
        else %if not, add mouse to end
            Mouse = Mouse(ind);
            Folder = Folder(ind);
            A1 = repmat('empty', 1,1);
            A2 = repmat('empty', 1,1);
            A3 = repmat('empty', 1,1);
            AddMouse = {Mouse, Folder, A1, A2, A3};

            RecordingOverview = [RecordingOverview; AddMouse];
         end
    end
end
    
%%
%start movement checking
for indMouse = 1:size(Mice, 1) %go per mouse
    Mouse = Mice(indMouse).name;
    
    for Acquisition = 1:3 %go per acquisition
        Recordings = {};
        
        for ind = 1:size(AcqList,2) %check per folder
            
            if  isequal(RecordingOverview{indMouse, (Acquisition+2)}, {'empty'}) && ... %check if not already done
                isequal(AcqList(ind).name(end-8:end-6), Mouse) && ... %check if it's the right mouse
                isequal(AcqList(ind).name(end-3:end-3), num2str(Acquisition)) && ... %check if it's the right acquisition
                isequal(AcqList(ind).name(end), '1') %get R1, see which acquisitions match
                
                %start to search all recordings of one acquisition, one mouse
                ToSearch = AcqList(ind).name(1:end-2);
                DoubleAcq = strfind(struct2cell(AcqList), ToSearch);
                DoubleAcq = reshape(DoubleAcq, 1, []);
                
                DoubleAcq(cellfun(@(x) isempty(x),DoubleAcq))= {0};
                DoubleAcq = cell2mat(DoubleAcq);
                
                for index = 1:size(AcqList,2)
                    if DoubleAcq(index) == 1 %if the file is marked as doubleacq or the R1
                        Recordings{end+1} = AcqList(index).name;
                        
                        findfileseps = strfind(AcqList(ind).name, filesep);
                        savefolder = AcqList(ind).name(1:findfileseps(end));
                        savename = AcqList(index).name(findfileseps(end-1): end);
                        savename = replace(savename, filesep, '-');
                        
                        if exist([savefolder 'Movement' savename '.tif'], 'file')
                            % if we already saved it, just show that one
                            MovIm = imread([savefolder 'Movement' savename '.tif']);
                            figure('Position', [50 200 400 400]);
                            imagesc(MovIm)
                        else
                            % otherwise, plot the movement
                            PlotTreadmill(AcqList(index).name) %plot R1, R2, R3 etc.
                            saveas(gcf, [savefolder 'Movement' savename], 'tiff');                 
                        end
                    end
                end

                [indx, ~] = listdlg('PromptString', {'Which recording do you want to use'},...
                    'SelectionMode','single','ListString', Recordings, 'ListSize', [600 300]);
                
                % add right recording to the overview
                row = matches(RecordingOverview.Mouse, Mouse);
                RecordingOverview(row, (Acquisition+2)) = Recordings(indx);
                close all
                
            end
            Recordings = {}; %to start for new acquisition
        end
    end
    save( '/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat', 'RecordingOverview')
end
end



















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
