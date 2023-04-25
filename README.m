%% README

% To start the pipeline, we want to have an overview of which mice and
% which acquisitions we have. We make this with MakeAcqList. The path in
% this script is hardcoded!!! Change if you have a different path.

[Mice, AcqList] = MakeAcqList;

% We then do the Pipeline for each mouse individually. This is the
% preprocessing, filtering etc., where the mice are not compared to each
% other yet. For more detailed descriptions, see inside the pipeline code.

for ind =1:size(AcqList,2)
    Pipeline_CaCl_GCaMP(AcqList(ind).name,1)
end