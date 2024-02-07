% Do first:
% Preprocessing_Pipeline_CaCl_GCaMP
% Umit_Pipeline_CaCl_GCaMP

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
        CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group'});
%         CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
%         
        %% Seed pixel correlation map
%         CombinedSPCM(datanames{ind1}, Acquisitions{ind2}, {'Group'}, 1, 1); %GSR = 1
        
        %% seed spread
%         CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
%         CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group'});
        
        %% GCaMP fluctuations
        CombinedGCaMPFluctuations(Acquisitions{ind2}, datanames{ind1}, {'Group', 'Sex'},Overwrite);
        CombinedGCaMPFluctuations_dotplot(Acquisitions{ind2}, datanames{ind1});
%         
        %% Connectivity percentage
    end
end

%% Movement plot
[pvalues] = CombinedMovementPlot({'Combi'}, 'Both', 0);
disp(pvalues.GCaMP)
disp(pvalues.Speckle)


