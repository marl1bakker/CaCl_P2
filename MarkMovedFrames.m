%     %% step 2; Mark frames with movement
%     aiFilesList = dir([DataFolder 'ai_*.bin']);
%     
%     % Opening of the files:
%     AnalogIN = [];
%     for ind = 1:size(aiFilesList,1)
%         data = memmapfile([DataFolder aiFilesList(ind).name],...
%             'Offset', 5*4, 'Format', 'double', 'repeat', inf);
%         tmp = data.Data;
%         tmp = reshape(tmp, 1e4, 11, []);
%         tmp = permute(tmp,[1 3 2]);
%         tmp = reshape(tmp,[],11);
%         AnalogIN = [AnalogIN; tmp];
%     end
%     
%     Treadmill = AnalogIN(:,5);
%     
%     % plot(Treadmill);
%     % line([1 6000000], [1.64, 1.64], 'Color', '#A2142F');
%     % line([1 6000000], [1.6455, 1.6455], 'Color', '#A2142F');
%     
%     Treadmill = (Treadmill < 1.6455) & (Treadmill > 1.64); % get a logical array for movement (0) or not (1)
%     
%     load([DataFolder 'green.mat']);
%     TimestepsPerFrame = size(Treadmill,1)/datLength;
%     
%     FramesMoved = ones(datLength,1);
%     ind1 = 1; %just to get it started
%     for ind = 1:datLength
%         ind2 = ceil(ind * TimestepsPerFrame);
%         
%         if any(Treadmill(ind1:ind2)<1) % gives 1 if there was movement
%             FramesMoved(ind) = 0; %Give a 0 for frames that have movement
%         end
%         
%         ind1 = floor(ind * TimestepsPerFrame) ; %this is already updated for next step
%     end
%     
%     clear aiFilesList AnalogIN Datatype datFile datLength datName datSize ...
%         dim_names FirstDim Freq ind ind1 ind2 Stim tExposure TimestepsPerFrame ...
%         tmp Treadmill