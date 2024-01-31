%% bug patch when mouse multiplied
%load .mat file of project
% load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\M30-35.mat')
% load('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\project_17-1-24.mat')

alreadydone = 'x';
newdirectoryPairs = [];

for ind = 1:size(directoryPairs, 2)
    name = directoryPairs{ind}{1};
    if sum(contains(alreadydone, name)) == 0
        newdirectoryPairs = [newdirectoryPairs, directoryPairs(ind)];
        alreadydone = [alreadydone; {name}];
    else
        disp('Mouse double')
    end
end

directoryPairs = newdirectoryPairs;

%save .mat file of project
% save('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\M30-35.mat', 'directoryPairs');
% save('C:\Users\marle\OneDrive\Documenten\MATLAB\P2\USGui\project_17-1-24.mat', 'directoryPairs');

%% check why it doesnt want to save

for ind = 1:size(directoryPairs,2)
    currentdata = directoryPairs{1,ind}{1,4};
    currentmouse = directoryPairs{1,ind}{1,1};

    if currentdata.done ~= 1 || size(currentdata.mean_diameter, 1) < 1 ...
            || size(currentdata.mean_velocity, 1) < 1
        disp(['stuck at mouse ' currentmouse])
    else
        try
            (max(currentdata.mean_diameter)-min(currentdata.mean_diameter));
            (max(currentdata.mean_velocity)-min(currentdata.mean_velocity))/mean(currentdata.mean_velocity);
            (max(currentdata.mean_velocity)-min(currentdata.mean_velocity))/max(currentdata.mean_velocity);
            (max(currentdata.mean_diameter)-min(currentdata.mean_diameter))/mean(currentdata.mean_diameter);
            mean(currentdata.mean_velocity);

        catch
            disp(['stuck at mouse ' currentmouse num2str(ind)])

        end
    end

end






