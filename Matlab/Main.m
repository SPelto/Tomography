% Create a sinogram
clear
close all
visualize = false;

bonePath = '/home/grozz/DATA';
% bonePath = '/Users/anttikoivurova/Documents/MATLAB/tomoproject/Koivurova_Peltonen/20191113_bone/corrected'; % mac
% bonePath = 'D:\gitProjects\Tomography\Data\20191113_bone\corrected'; % windows
boneName = '20191113_bone';

% Wanted sinogram parameters
binning = 4;
shift = 2;
N = 2048 / binning;

% We want the filename to have information about the angles, hence we give
% the angle as a string so that it can be easily used in the filename.
angleString = '90:6:269'; % Now with normal angle indexing

% Add angles as a variable based on the string
eval(['angles = ' angleString ';']) 

% Filenames doesn't allow for ":" characters, so we'll change them into "_"
angleString(angleString == ':') = '_';

sinoFileName = ['Sinogram_' angleString '_bin_'  num2str(binning) '_shift_' num2str(shift) '.mat'];

%% Check if we have sinogram with wanted parameters already
if ~isfile(sinoFileName)
    sinogram = createSino(bonePath, boneName, binning, shift, angles, visualize);
    figure('name', 'Sinogram')
    imshow(sinogram, [])
    save(sinoFileName, 'sinogram')
else 
    load(sinoFileName)
end




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
recn = reshape(A.'*sinogram_T(:),[N,N]);

% Plot the reconstruction
% figure(66)
% clf
% imagesc(recn)
% colormap gray

%% Reconstruction with Tikhonov
close all
alphas = [20:20:100];
X = zeros(1,length(alphas));
Y = zeros(1,length(alphas));

for iii=1:length(alphas)
    alpha = alphas(iii);
    Tik_rec = Solve_Tikh(sinogram_T, A, N, alpha);
    figure;
    imshow(Tik_rec, [], 'InitialMagnification', 200);
    Tik_rec_normalised = Tik_rec(:) / max(Tik_rec(:));
    
    Tik_rec_filename = ['Tik_rec_sparse_alpha_' alpha '_angles_' angleString '_binning_' num2str(binning) '_shift_' num2str(shift)];
    save(Tik_rec_filename, 'Tik_rec_normalised');
    
    X(iii) = log(norm(A*Tik_rec(:) - sinogram_T(:)));
    Y(iii) = log(norm(Tik_rec(:)));
end

figure;
plot(X,Y, 'bo-');
Xn = (X- min(X)) / (max(X)-min(X));
Yn = (Y- min(Y)) / (max(Y)-min(Y));
[~, Tik_min] = min(Xn.^2 + Yn.^2);


%% Total variation

alphas = [0.01];
X = zeros(1,length(alphas));
Y = zeros(1,length(alphas));

for iii=1:length(alphas)
    alpha = alphas(iii);
    TV_rec_sparse = Solve_TV2(sinogram_T, A, N, alpha);
    figure;
    imshow(TV_rec_sparse, [], 'InitialMagnification', 200);
    title("alpha: " + alpha)
end

TV_rec_filename = ['TV_rec_sparse_alpha_001' '_angles_' angleString '_binning_' num2str(binning) '_shift_' num2str(shift)];
save(TV_rec_filename, 'TV_rec_sparse');
% 
% 
% %% Calculate reconstruction using ifanbeam
% DistanceSourceDetector  = 552.18;
% DistanceOffsetSample    = 275; % Bone
% % DistanceOffsetSample    = 75;  % Tomato
% DistanceSourceOrigin    = 109.83 + DistanceOffsetSample;
% pixelSize               = 0.05*binning;
% M                       = DistanceSourceDetector / DistanceSourceOrigin;
% effPixelSize            = pixelSize / M;
% D                       = DistanceSourceOrigin / effPixelSize;
% 
% 
% filter = ["Ram-Lak"];
% % filter = ["Shepp-Logan"];
% % filter = ["Cosine"];
% % filter = ["Hamming"];
% % filter = ["Hann"];
% % filter = ["None"];
% 
% reconstruction = ifanbeam(sinogram, D,...
%     'FanSensorGeometry', 'line',...
%     'OutputSize', 2240 / binning,...
%     'Filter', filter);
% 
% figure(1)
% clf
% imshow(reconstruction,[])
% title(filter)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [-1, 0, 1, 1]); % Sets the figure to open at the left screen.




