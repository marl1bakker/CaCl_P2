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
ylim([1.63 1.65])
title(DataFolder(end-9:end-1))
%to set lines, but the range changes per acq
% line([1 size(Treadmill,1)], [UpperLimit UpperLimit], 'Color', 'red') 
% line([1 size(Treadmill,1)], [LowerLimit LowerLimit], 'Color', 'red')

NotMoved = sum((Treadmill < UpperLimit) & (Treadmill > LowerLimit)); % get a logical array for movement (0) or not (1)
Moved = size(AnalogIN,1) - NotMoved;
FractionMoved = Moved/size(AnalogIN,1);

subtitle(['Fraction movement = ' num2str(FractionMoved)])

end