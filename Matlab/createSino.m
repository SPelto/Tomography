function sinogram = createSino(filePath, photoName, binning, shift, visualize)
    addpath(filePath);

    bkgSize = 2^7;

    sinogram = zeros(2240,360);

    figure(1)
    clf
    for k = 1:360
        if mod(k,10) == 0
            disp(k)
        end
        fileName = [photoName '_' sprintf('%03d', k) '.tif'];
        im = double(imread(fileName));
        background = im(1:bkgSize/binning,1:bkgSize/binning);

        I0 = mean(background(:));
        I = im(end/2,:);

        lineIntegral = -log(I./I0);
        sinogram(:,k) = lineIntegral;

        if visualize
            sinogram2 = sinogram;
            sinogram2(:,k) = max(sinogram(:));
            
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

