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

% We want the filename to have information about the angles, hence we give
% the angle as a string so that it can be easily used in the filename.
angleString = '2:2:361'; % Not 1:360 because the first picture measured is a blank one

% Add angles as a variable based on the string
eval(['angles = ' angleString ';']) 

% Filenames doesn't allow for ":" characters, so we'll change them into "_"
angleString(angleString == ':') = '_';

sinoFileName = ['Sinogram_' angleString '_bin_'  num2str(binning) '_shift_' num2str(shift) '.mat'];
% Check if we have sinogram with wanted parameters already
if ~isfile(sinoFileName)
    sinogram = createSino(bonePath, boneName, binning, shift, angles, visualize);
    figure('name', 'Sinogram')
    imshow(sinogram, [])
    save(sinoFileName, 'sinogram')
else 
    load(sinoFileName)
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
A_filename = ['A_angles_' angleString '_binning_' num2str(binning) '_shift_' num2str(shift)];

if ~isfile(A_filename)
    A = createSysMat(binning,angles);
    save(A_filename, 'A')
else
    load(A_filename)
end
% ASTRA uses different orientation for the sinogram
sinogram_T = sinogram.';
recn = reshape(A.'*sinogram_T(:),[2048/binning,2048/binning]);

% Plot the reconstruction
figure(66)
clf
imagesc(recn)
colormap gray










