%% Upload to DANDI - NWB
% Get data here: https://doi.org/10.48324/dandi.001210/0.241111.1757

% to transform .dat files to .nwb files and vice versa
% direction can be nwb-dat or dat-nwb or bin-nwb (calls
% ImagesClassification)
% MouseTable only necessary for dat/bin-nwb, not the other way around
% load('/home/mbakker/P2_scripts/MarleenP2/MiceCodes.mat')

% best would be to do it from .bin files directly like:
% DANDI_NWB_conversion('/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2', 'bin-nwb', '/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2/bin-nwb', Mice)

% DANDI_NWB_conversion('/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2', 'dat-nwb', '/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2/dat-nwb')
% DANDI_NWB_conversion('/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2/dat-nwb', 'nwb-dat', '/home/mbakker/P2_scripts/MarleenP2/ExampleData/M32-A1-R2/nwb-dat')

function DANDI_NWB_conversion(DataFolder, direction, SaveFolder, MouseTable)

if( ~strcmp(DataFolder(end), filesep) )
    DataFolder = [DataFolder filesep];
end

if ~exist('SaveFolder', 'var')
    SaveFolder = DataFolder;
elseif( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = [SaveFolder filesep];
end

if ~exist('direction', 'var')
    direction = 'dat-nwb';
end

%% set up
switch direction
    case {'dat-nwb', 'bin-nwb'}
        datanames = {'fluo_567','red','green'};
        descriptions = {'GCaMP', 'Red', 'Green'};
        emissions = [567., 620., 535.];
        excitations = [472., 620., 535.];

        % Get mouse and acquisition etc. PARTIALLY HARDCODED
        Mouse_ID = DataFolder(regexp(DataFolder, 'M[123]\w*-A[123]-R\w*'):end);
        indseps = strfind(Mouse_ID, filesep);
        Mouse_ID = Mouse_ID(1:indseps(1)-1);
        % indseps = strfind(DataFolder, filesep);
        % Mouse_ID = DataFolder(indseps(end-1)+1:indseps(end)-1);
        inddash = strfind(Mouse_ID, '-');
        Mouse = Mouse_ID(1:inddash(1)-1);
        if exist('MouseTable', 'var')
            indmouse = find(matches(MouseTable.CodeOfMouse, Mouse));
            sex_mouse = char(MouseTable.MaleFemale(indmouse));
            sex_mouse = sex_mouse(1);
            mouse_group = char(MouseTable.CaClSham(indmouse));
            age_mouse = 'P14W/P19W'; %14 to 19 weeks ish
        else
            sex_mouse = 'U';
            mouse_group = 'Unknown';
            age_mouse = 'U';
        end
        clear indseps inddash indmouse MouseTable

        if matches(direction, 'bin-nwb')
            ImagesClassification(DataFolder, SaveFolder, 1, 1, 0);
            MarkMovedFrames(DataFolder, SaveFolder);
            DataFolder = SaveFolder;
        end

        % get general information
        load([DataFolder 'AcqInfos.mat'], 'AcqInfoStream');
        dt = AcqInfoStream.DateTime;
        dt(strfind(dt, '_')) = [];
        dt = datetime(str2double(dt(1:4)), str2double(dt(5:6)), str2double(dt(7:8)), str2double(dt(9:10)), str2double(dt(11:12)), str2double(dt(13:14)), 'TimeZone', 'local'); %prob better ways to do it

        % subject information
        subject = types.core.Subject(...,
            'subject_id', Mouse,...
            'description', mouse_group,...
            'sex', sex_mouse,...
            'age', age_mouse,...
            'species', 'Mus musculus');

        % device information
        device = types.core.Device(...
            'description', ['Camera Model: ' AcqInfoStream.Camera_Model ' on LightTrackOiS200 system'],...
            'manufacturer', 'Labeo Technologies');

    case 'nwb-dat'
        datanames = dir(DataFolder);
        datanames = struct2table(datanames);
        datanames = datanames(~datanames.isdir,:).name;
        for indnames = 1:size(datanames,1)
            if ~strcmp(datanames{indnames}(end-3:end), '.nwb')
                datanames(indnames) = {'x'};
            end
        end
        datanames = datanames(~matches(datanames, 'x'));
end

%% go per channel
for ind = 1:3

    datname = datanames{ind};

    switch direction
        % save as .nwb file
        case {'dat-nwb', 'bin-nwb'}

            load([DataFolder datname '.mat'], 'Freq', 'datSize', 'tExposure')

            % general nwb file
            nwb = NwbFile( ...
                'session_description', 'combined GCaMP IOI data',...
                'identifier', [Mouse_ID '-' descriptions{ind}], ...
                'session_start_time', dt,...
                'general_experimenter', 'Bakker, Marleen', ...
                'general_institution', 'Institut Cardiologie Montreal');
            % 'general_related_publications', {'DOI:...'}); % optional

            % subject
            nwb.general_subject = subject;

            % device
            nwb.general_devices.set('Device', device);

            % imaging plane
            optical_channel = types.core.OpticalChannel(...
                'description', descriptions(ind), ...
                'emission_lambda', emissions(ind));
            imaging_plane_name = 'imaging_plane';
            imaging_plane = types.core.ImagingPlane( ...
                'optical_channel', optical_channel, ...
                'description', 'whole brain', ...
                'device', types.untyped.SoftLink(device), ...
                'excitation_lambda', excitations(ind), ...
                'imaging_rate', Freq, ...
                'indicator', descriptions(ind), ...
                'location', 'cortex');

            nwb.general_optophysiology.set(imaging_plane_name, imaging_plane);

            % get .dat data
            fid = fopen([DataFolder datname '.dat']);
            dat = fread(fid, inf, '*single');
            dat = reshape(dat, datSize(1), datSize(2), []);
            fclose(fid);

            % first dimension should be time in .nwb files:
            dat = permute(dat, [3 1 2]);
            %imagesc(reshape(dat(1,:,:), 512, 512));

            % store as one photon data
            InternalOnePhoton = types.core.OnePhotonSeries( ...
                'data', dat, ...
                'description', datname,...
                'imaging_plane', types.untyped.SoftLink(imaging_plane), ...
                'starting_time', 0., ...
                'starting_time_rate', Freq,...
                'binning', AcqInfoStream.Binning,... % is this same?
                'exposure_time', tExposure,...
                'data_unit', '-');
            nwb.acquisition.set('1pInternal', InternalOnePhoton);

            % store the movement data and acq info with first data you save:
            if ind == 1
                load([DataFolder 'MovMask.mat'], 'MovMask')

                spatial_series_ts = types.core.SpatialSeries(...
                    'data', MovMask, ...
                    'starting_time', 0., ...
                    'starting_time_rate', Freq...
                    );
                Position = types.core.Position('SpatialSeries', spatial_series_ts);
                behavior_mod = types.core.ProcessingModule('description', 'Additional Acquisition information and movement measured by treadmill, already in binary (0 is movement, 1 is no movement)');
                behavior_mod.nwbdatainterface.set('Position', Position);
                nwb.processing.set('behavior', behavior_mod);
            end

            % export
            nwbExport(nwb, [SaveFolder datname '.nwb']);
            clear nwb optical_channel dat fid InternalOnePhoton





        case 'nwb-dat'
            % open nwb file
            if ~strcmp(datname(end-3:end), '.nwb')
                datname = [datname '.nwb'];
            end
            readnwb = nwbRead([DataFolder datname], 'ignorecache');
            savedatname = readnwb.acquisition.get('1pInternal').description;

            % check if there is a movement mask in this file
            try 
                MovMask = readnwb.processing.get('behavior').nwbdatainterface.get('Position').spatialseries.get('SpatialSeries').data.load;
                MovMask = MovMask';
                save([SaveFolder 'MovMask.mat'], 'MovMask')
                % save([SaveFolder 'AcqInfos.mat'], 'AcqInfoStream')
            catch
                % disp('No movmask in this file.')
            end

            % Get data
            dat = readnwb.acquisition.get('1pInternal').data.load;
            dat = permute(dat, [2, 3, 1]); % turn back into time last

            % save as dat
            fid = fopen([SaveFolder savedatname '.dat'],'w');
            fwrite(fid,dat,'single');
            fclose(fid);

            % make corresponding .mat file
            Datatype = 'single';
            FirstDim = 'y'; %hardcoded
            Freq = readnwb.general_optophysiology.get('imaging_plane').imaging_rate;
            datFile = [savedatname '.dat'];
            datLength = length(dat);
            datName = 'data'; %hardcoded
            datSize = [size(dat,1), size(dat,2)];
            dim_names = {'Y','X','T'}; %hardcoded
            tExposure = readnwb.acquisition.get('1pInternal').exposure_time;

            save([SaveFolder savedatname '.mat'], 'Datatype','FirstDim', 'Freq', ...
                'datFile', 'datLength', 'datName', 'datSize', 'dim_names', 'tExposure')

            clear dat readnwb fid MovMask
    end
end

%% make AcqInfos.mat & fluo 475.dat
% this does not save as much in the acqinfos.mat file as the original
% pipeline. 
switch direction
    case 'nwb-dat'
        readnwb = nwbRead([DataFolder datname], 'ignorecache'); % does not matter which you pick, so pick last one
        load([SaveFolder savedatname '.mat'], 'datSize')

        AcqInfoStream.DateTime = char(readnwb.session_start_time(1));
        AcqInfoStream.FrameRateHz = 30; %hardcoded
        AcqInfoStream.Width = datSize(1);
        AcqInfoStream.Height = datSize(2);
        AcqInfoStream.Binning = readnwb.acquisition.get('1pInternal').binning;

        cameradescription = readnwb.general_devices.get('Device').description;
        camera_id_string = 'Camera Model: '; %what you put in the description before camera model
        indstart = strfind(cameradescription, camera_id_string);
        indstart = indstart(1)+length(camera_id_string);
        indend = strfind(cameradescription(indstart:end), ' ');
        indend = indend(1)+length(camera_id_string)-1;
        AcqInfoStream.Camera_Model = cameradescription(indstart:indend);
        clear indstart indend camera_id_string cameradescription
        AcqInfoStream.Stimulation = 0;
        AcqInfoStream.MultiCam = 1;

        % Illumination channels HARDCODED
        % ID Color CamIdx FrameIdx
        AcqInfoStream.Illumination1 = struct('ID', 1, 'Color', 'Red', 'CamIdx', 2, 'FrameIdx', 1);
        AcqInfoStream.Illumination2 = struct('ID', 2, 'Color', 'Green', 'CamIdx', 1, 'FrameIdx', 1);
        AcqInfoStream.Illumination3 = struct('ID', 3, 'Color', 'Fluo #1 475 nm', 'CamIdx', 2, 'FrameIdx', 2);
        AcqInfoStream.Illumination4 = struct('ID', 4, 'Color', 'Fluo #2 567 nm', 'CamIdx', 1, 'FrameIdx', 2);

        save([SaveFolder 'AcqInfos.mat'], 'AcqInfoStream');

        % make a fluo 475.dat and .mat to not fuck up the correctdichroic script
        fid = fopen([SaveFolder 'fluo_475.dat'],'w');
        fwrite(fid,zeros(datSize, datLength),'single'); %take size from prev.
        fclose(fid);

        Datatype = 'single';
        FirstDim = 'y'; %hardcoded
        Freq = 15; %hardcoded
        datFile =  'fluo_475.dat';
        datLength = length(dat);
        datName = 'data'; %hardcoded
        datSize = [512, 512];
        dim_names = {'Y','X','T'}; %hardcoded
        tExposure = readnwb.acquisition.get('1pInternal').exposure_time;

        save([SaveFolder 'fluo_475.mat'], 'Datatype','FirstDim', 'Freq', ...
            'datFile', 'datLength', 'datName', 'datSize', 'dim_names', 'tExposure')




        % delete .dat files if you go bin - nwb directly
    case 'bin-nwb'
        del_answer = questdlg('Delete .dat files?');
        if matches(del_answer, 'Yes')
            delete([SaveFolder 'AcqInfos.mat'])
            delete([SaveFolder 'fluo_475.dat'])
            delete([SaveFolder 'fluo_475.mat'])
            delete([SaveFolder 'MovMask.mat'])
            for ind = 1:3
                delete([SaveFolder datanames{ind} '.dat'])
                delete([SaveFolder datanames{ind} '.mat'])
            end
        end
end
end
