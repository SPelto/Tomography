function sinogram = createSino(filePath, photoName, binning, shift, angles, visualize)
    addpath(filePath);
    
    bkgSize = 2^8;

    sinogram = zeros(2240 / binning, length(angles));

    figure(1)
    clf
    
    ind = 0;
    for k = angles
        ind = ind + 1;
        % Display progress
        if mod(k,10) == 0
            disp([k length(angles)])
        end
        
        % Load the picture
        fileName = [photoName '_' sprintf('%03d', k) '.tif'];
        im = double(imread(fileName));
        
        % Shift the picture to center it based on rotation axis
        im = circshift(im, shift, 1);

        % Do binning
        if binning > 1
            im = bin_projection(im, binning);
        end
        background = im(1:bkgSize/binning,1:bkgSize/binning);
        I0 = mean(background(:));
        I = im(end/2,:);

        lineIntegral = -log(I./I0);
        
        sinogram(:,ind) = lineIntegral;

        if visualize
            sinogram2 = sinogram;
            sinogram2(:,k - 1) = max(sinogram(:));
            
            im2 = im;
            im2(end/2,:) = max(im(:));
            subplot(1,2,1)
            imshow(im2,[])
            colormap gray

            subplot(1,2,2)
            imagesc(sinogram2)
            colormap gray
            pause(0.1)
        end
    end
end

