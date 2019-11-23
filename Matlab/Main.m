% Create a sinogram
clear
close all
visualize = false;

bonePath = '/Users/anttikoivurova/Documents/MATLAB/tomoproject/Koivurova_Peltonen/20191113_bone/corrected'; % mac
% bonePath = 'D:\Tomography Project\Data\20191113_bone\corrected'; % windows
boneName = '20191113_bone';

tomatoPath = '/Users/anttikoivurova/Documents/MATLAB/tomoproject/Koivurova_Peltonen/20191113_tomato/corrected'; % mac
% tomatoPath = 'D:\Tomography Project\Data\20191113_tomato\corrected'; % windows 
tomatoName = '20191113_tomato';

% Wanted sinogram parameters
binning = 1;
shift = 0;
angles = 2:361; % Not 1:360 because the first picture measured is a blank one

% Check if we have sinogram with wanted parameters already
if ~isfile("Sinogram_full_bin_" + num2str(binning) + "_shift_" + num2str(shift) + ".mat")
    sinogram = createSino(bonePath, boneName, binning, shift, angles, visualize);
    figure('name', 'Sinogram')
    imshow(sinogram, [])
    sinoFileName = "Sinogram_full_bin_" + num2str(binning) + "_shift_" + num2str(shift);
    save(sinoFileName, 'sinogram')
else 
    load("Sinogram_full_bin_" + num2str(binning) + "_shift_" + num2str(shift))
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




% filters = ["Ram-Lak" "Shepp-Logan" "Cosine" "Hamming" "Hann" "None"]
% Choose a single filter if you do not wish to iterate over all
% possibilities

% filters = ["Ram-Lak"];
% filters = ["Shepp-Logan"];
% filters = ["Cosine"];
% filters = ["Hamming"];
% filters = ["Hann"];
filters = ["None"];

ind = 0;
for filter = filters
    ind = ind + 1;
    reconstruction = ifanbeam(sinogram, D,'FanSensorGeometry', 'line', 'OutputSize', 2240 / binning, 'Filter', filter);
    figure(ind)
    clf
    imshow(reconstruction,[])
    title(filter)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [-1, 0, 1, 1]);
end






