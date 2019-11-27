% Create a sinogram
clear
close all
visualize = false;

% bonePath = '/Users/anttikoivurova/Documents/MATLAB/tomoproject/Koivurova_Peltonen/20191113_bone/corrected'; % mac
bonePath = 'D:\gitProjects\Tomography\Data\20191113_bone\corrected'; % windows
boneName = '20191113_bone';

% tomatoPath = '/Users/anttikoivurova/Documents/MATLAB/tomoproject/Koivurova_Peltonen/20191113_tomato/corrected'; % mac
% tomatoPath = 'D:\gitProjects\Tomography\Data\20191113_tomato\corrected'; % windows 
tomatoName = '20191113_tomato';

% Wanted sinogram parameters
binning = 4;
shift = 0;
angles = 2:1:361; % Not 1:360 because the first picture measured is a blank one

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
pixelSize               = 0.05*binning;
M                       = DistanceSourceDetector / DistanceSourceOrigin;
effPixelSize            = pixelSize / M;
D                       = DistanceSourceOrigin / effPixelSize;


filter = ["Ram-Lak"];
% filter = ["Shepp-Logan"];
% filter = ["Cosine"];
% filter = ["Hamming"];
% filter = ["Hann"];
% filter = ["None"];

reconstruction = ifanbeam(sinogram, D,...
    'FanSensorGeometry', 'line',...
    'OutputSize', 2240 / binning,...
    'Filter', filter);

figure(1)
clf
imshow(reconstruction,[])
title(filter)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [-1, 0, 1, 1]); % Sets the figure to open at the left screen.


%% Create the forward projector
A = createSysMat(binning);
whos A

% ASTRA uses different orientation for the sinogram
sinogram_T = sinogram.';

recn = reshape(A.'*sinogram_T(:),[2048/binning,2048/binning]);

% Plot the reconstruction
figure(66)
clf
imagesc(recn)
colormap gray
















