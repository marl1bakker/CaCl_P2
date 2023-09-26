function [Done] = CheckIfDone(SaveFolder, FunctionToCheck)

if( ~strcmp(SaveFolder(end), filesep) )
    SaveFolder = [SaveFolder filesep];
end

load('/media/mbakker/GDrive/P2/GCaMP/LogBook.mat'); %hardcoded
Mouse = SaveFolder(end-16:end-14);
Acq = SaveFolder(end-12:end-8);

try
    inds = find(strcmp(Mouse, LogBook.Subject));
    LogBook = LogBook(inds,:);
    inds = find(strcmp(Acq, LogBook.Acquisition));
    LogBook = LogBook(inds,:);
    inds = find(strcmp(FunctionToCheck, LogBook.TaskName));
    LogBook = LogBook(inds,:);
    inds = find(1, LogBook.Completed);
    LogBook = LogBook(inds,:);
catch
    %if it errors any of the above checks, it's good to go
    disp('Check function CheckIfDone, might have something weird')
end

if size(LogBook,1) > 0 %to do last check, see if it's already done
    disp([FunctionToCheck ' ' Mouse ' ' Acq ' already run'])
    Done = 1;
else 
    Done = 0;
end
end