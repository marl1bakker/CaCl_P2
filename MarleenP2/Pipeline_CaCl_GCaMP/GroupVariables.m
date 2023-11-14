% Give the grouping variables in Grouping. This can be Sex, Group, or
% {'Group', 'Sex'};. The function will add a column to the
% RecordingOverview with your Combi if your grouping exists of more than
% one. It will also return the variable groups, which is a name that is an
% addition of all group variables, and in the second column a code unique
% for that group. With this code, you can prevent double makings of plots
% etc. For now, the function will give either 1 or 2 per group variable.
% This will not work with more than 2 groups for now. 


function [groups, RecordingOverview] = GroupVariables(RecordingOverview, Grouping)

%% Make a group column in RecordingOverview or use an existing one. 
switch size(Grouping,2)
  
    case 0
        disp('No grouping variable selected, function exited')
        return
        
    case 1
        eval(['currentgroup = cellstr(RecordingOverview.' Grouping{1} ');' ]);
        RecordingOverview.Codes = cellstr(repmat(' ', size(RecordingOverview,1),1));
        
        %codes
        codes = zeros(size(RecordingOverview,1),1);
        groupoptions = string(unique(currentgroup));
        for ind2 = 1:size(groupoptions, 1)
            codes = codes + (ind2 .* contains(currentgroup, groupoptions{ind2}));
        end
        RecordingOverview.Codes = strcat(RecordingOverview.Codes, num2str(codes));
        
        %groups
        groups = [cellstr(unique(currentgroup)) unique(RecordingOverview.Codes)];
        RecordingOverview.Codes = str2double(RecordingOverview.Codes);
        
    otherwise %more than 1 grouping variable - make a Combi group
        RecordingOverview.Combi = cellstr(repmat(' ', size(RecordingOverview,1),1));
        RecordingOverview.Codes = cellstr(repmat(' ', size(RecordingOverview,1),1));
        codes = zeros(size(RecordingOverview,1),1);
        
        for ind = 1:size(Grouping, 2)
            %full names
            eval(['currentgroup = cellstr(RecordingOverview.' Grouping{ind} ');' ]);
            RecordingOverview.Combi = strcat(RecordingOverview.Combi, currentgroup);
            
            %codes
            groupoptions = string(unique(currentgroup));
            for ind2 = 1:size(groupoptions, 1)
                codes = codes + (ind2 .* contains(currentgroup, groupoptions{ind2}));
            end
            RecordingOverview.Codes = strcat(RecordingOverview.Codes, num2str(codes));
            codes = zeros(32,1);
        end
        
        groups = [unique(RecordingOverview.Combi) unique(RecordingOverview.Codes)];
%         Grouping = 'Combi';
        RecordingOverview.Combi = categorical(RecordingOverview.Combi);
        RecordingOverview.Codes = str2double(RecordingOverview.Codes);
end

end