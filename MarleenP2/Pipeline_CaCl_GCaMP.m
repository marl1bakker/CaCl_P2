% Do first:
% Preprocessing_Pipeline_CaCl_GCaMP
% Umit_Pipeline_CaCl_GCaMP

%% Mice combined
% datanames = {'hemoCorr_fluo'};
% datanames = {'hemoCorr_fluo', 'HbO', 'HbR'};
% Acquisitions = {'A1','A2','A3'};
% % Acquisitions = {'A1'};
% Overwrite = 1; %don't overwrite stuff
% 
% for ind1 = 1:size(datanames,2)
%     for ind2 = 1:size(Acquisitions,2)
% 
%         disp(['***** ' datanames{ind1} ' ' Acquisitions{ind2} ' *****'])
% 
%         %% Correlation matrix
%         % CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group'});
% %         CombinedCorrMatrix(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
% %         
%         %% Seed pixel correlation map
%         % CombinedSPCM(datanames{ind1}, Acquisitions{ind2}, {'Group'}, 1, 1); %GSR = 1
% 
%         %% seed spread
% %         CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group', 'Sex'});
%         CombinedSeedSpread(datanames{ind1}, Acquisitions{ind2}, {'Group'});
% 
%         %% GCaMP fluctuations
%         CombinedGCaMPFluctuations(Acquisitions{ind2}, datanames{ind1}, {'Group', 'Sex'},Overwrite);
%         CombinedGCaMPFluctuations_dotplot(Acquisitions{ind2}, datanames{ind1});
% %         
%         %% Connectivity percentage
% 
%     end
% end


%% Figure 1, timecourse single pixel
ExampleData_GCaMP_HbO_HbR('/media/mbakker/GDrive2/P2/GCaMP/M32/A1-R2/CtxImg')

%% Figure 2&3 - NVC in sham
plotNVCinRS('A1', 'OutliersRemoved') %gives for both paper and appendix
load('/media/mbakker/GDrive/P2/GCaMP/NVC/Sham_RS/allSpecs.mat', 'allSpecs');
StatsNVCinRS(allSpecs)

%% Figure 4 - ultrasound (on Dell, not linux!!)
% load('C:\Users...)
% PlotUSResultsBoxplot(ResultsUS, 'Selection', 'R')

%% Figure 5 - Correlation Matrices
CorrMatricesAllAcq
CorrMatricesAllAcq('HbO')
CorrMatricesAllAcq('HbR')
[Stats] = CorrMatricesStats; %saves as excel file for appendix

%% Figure 6 - NVC in sham vs cacl
plotNVCAllAcq('R', 'OutliersRemoved') % in paper
plotNVCAllAcq('BigROI', 'OutliersRemoved') % in appendix


%% Movement plot (appendix)
[pvalues] = CombinedMovementPlot({'Combi'}, 'Both', 0);
disp(pvalues.GCaMP)
disp(pvalues.Speckle)

% [pvalues] = CombinedMovementPlot({'Group'}, 'Both', 0);
[pvalues] = CombinedMovementPlot({'Group'}, 'GCaMP', 0);
disp(pvalues.GCaMP)
% disp(pvalues.Speckle)









