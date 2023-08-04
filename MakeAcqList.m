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