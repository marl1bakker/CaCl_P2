% Do first:
% Preprocessing_Pipeline_CaCl_GCaMP
% Umit_Pipeline_CaCl_GCaMP

%% Set-up
% load('/home/mbakker/P2_scripts/MarleenP2/RecordingOverview.mat')
% Recordings = [RecordingOverview.A1; RecordingOverview.A2; RecordingOverview.A3];
% Mice = [RecordingOverview.Mouse; RecordingOverview.Mouse; RecordingOverview.Mouse];
% SaveDirs = [RecordingOverview.SaveDirectory; RecordingOverview.SaveDirectory;RecordingOverview.SaveDirectory;];
%
% Overwrite = 0; %don't overwrite stuff
%
% %% single mouse
% for ind = 1:size(Recordings,1)
%     Mouse = Mice{ind};
%     DataFolder = Recordings{ind};
%
%     if( ~strcmp(DataFolder(end), filesep) )
%         DataFolder = [DataFolder filesep];
%     end
%
%     SaveFolder = [SaveDirs{ind} filesep Mouse filesep DataFolder(end-5:end) 'CtxImg' filesep];
%     disp(Mouse)
%     disp(DataFolder(end-9:end-1))
%
%     %% Get Timecourses
%     GetTimecourses(SaveFolder, 'average', Overwrite, 'hemoCorr_fluo')
%     GetTimecourses(SaveFolder, 'average', Overwrite, 'HbO')
%     GetTimecourses(SaveFolder, 'average', Overwrite, 'HbR')
%
%     %% calculate the tform to align the mice/brains
%     AlignAllBrains(SaveFolder);
%
%     %% possible single subject figures
%     % Make corr Matrix single subject
%     %     SingleSubjectCorrMatrix(SaveFolder, Overwrite)
%
%     % Make SPCM single subject
%     %     SingleSubjectSPCM(SaveFolder, 'hemoCorr_fluo')
%
% end

%% Mice combined
% datanames = {'hemoCorr_fluo'};
datanames = {'hemoCorr_fluo', 'HbO', 'HbR'};
Acquisitions = {'A1','A2','A3'};
% Acquisitions = {'A1'};
Overwrite = 0; %don't overwrite stuff

for ind1 = 1:size(datanames,2)
    for ind2 = 1:size(Acquisitions,2)
        
        disp(['***** ' datanames{ind1} ' ' Acquisitions{ind2} ' *****'])
        
        %% Correlation matrix
        %         CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group'});
        %         CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
        
        %% Seed pixel correlation map
        %         CombinedSPCM(datanames{ind1}, Acquisitions{ind2}, {'Group'}, 1, 1); %GSR = 1
        %         CombinedSPCM(datanames{ind1}, Acquisitions{ind2}, 1, SaveDir, Overwrite); %GSR = 1
        
        %% seed spread
        CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
        
        % If you do the {'Group', 'Sex'} variables ones, you can do this code
        % after, to not have to calculate 'Groups':
        
        if ~exist(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} filesep  'CaCl_' Acquisitions{ind2} '.mat'], 'file')
            load(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} '/CaClFemale_' Acquisitions{ind2} '.mat']);
            load(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} '/CaClMale_' Acquisitions{ind2} '.mat']);
            spreadsCaCl = [spreadsCaClFemale, spreadsCaClMale];
            save(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} filesep  'CaCl_' Acquisitions{ind2} '.mat'], ...
                'spreadsCaCl'); %save the matrix
            clear spreadsCaCl*
        end
        
        if ~exist(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} filesep  'Sham_' Acquisitions{ind2} '.mat'], 'file');
            load(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} '/ShamFemale_' Acquisitions{ind2} '.mat']);
            load(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} '/ShamMale_' Acquisitions{ind2} '.mat']);
            spreadsSham = [spreadsShamFemale, spreadsShamMale];
            save(['/media/mbakker/GDrive/P2/GCaMP/SeedSpread/Combined/' datanames{ind1} filesep  'Sham_' Acquisitions{ind2} '.mat'], ...
                'spreadsSham'); %save the matrix
            clear spreadsSham*
        end
        
        CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group'});
        
    end
end
