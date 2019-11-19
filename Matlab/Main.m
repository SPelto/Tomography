% Create a sinogram
clear
visualize = true;
bonePath = 'D:\Tomography Project\Data\20191113_bone\corrected';
boneName = '20191113_bone';

tomatoPath = 'D:\Tomography Project\Data\20191113_tomato\corrected';
tomatoName = '20191113_tomato';

if ~isfile("Sinogram_full.mat")
    sinogram = createSino(bonePath, boneName, visualize);
    figure('name', 'Sinogram')
    imshow(sinogram, [])
    save('Sinogram_full', 'sinogram')
else 
    load('Sinogram_full')
end


%% Calculate reconstruction using ifanbeam
DistanceSourceDetector  = 552.18;
DistanceOffsetSample    = 275; % Bone
% DistanceOffsetSample    = 75;  % Tomato
DistanceSourceOrigin    = 109.83 + DistanceOffsetSample;
pixelSize               = 0.05;
M                       = DistanceSourceDetector / DistanceSourceOrigin;
effPixelSize            = pixelSize / M;
D                       = DistanceSourceOrigin / effPixelSize;




filters = ["Ram-Lak" "Shepp-Logan" "Cosine" "Hamming" "Hann" "None"]
% Choose a single filter if you do not wish to iterate over all
% possibilities

% filters = 'Ram-Lak'
% filters = 'Shepp-Logan'
% filters = 'Cosine'
% filters = 'Hamming'
% filters = 'Hann'
% filters = 'None'

ind = 0;
for filter = filters
    ind = ind + 1;
    reconstruction = ifanbeam(sinogram, D,'FanSensorGeometry', 'line', 'OutputSize', 2240, 'Filter', filter);
    figure(ind)
    clf
    imshow(reconstruction,[])
    title(filter)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [-1, 0, 1, 1]);
end






