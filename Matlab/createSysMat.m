function A = createSysMat(dims,binning) 
    
   %%%% Creates forward projector operator.
   % dims       = Size of the image
   % binning    = binningfactor
   
    
    DistanceOffsetSample    = 275; % Bone
    DistanceSourceDetector  = 552.18;
    DistanceSourceOrigin    = 109.83 + DistanceOffsetSample;
    DistanceOriginDetector  = DistanceSourceDetector - DistanceSourceOrigin;
    pixelSize               = 0.05;
    M                       = DistanceSourceDetector / DistanceSourceOrigin;
    effPixelSize            = pixelSize / M;
    D                       = DistanceSourceOrigin / effPixelSize;
    
    
    % I think the problem lies in selecting these.
    detector_geom   = 'fanflat';
    det_width       = 1%effPixelSize; % Especially this one
    det_count       = 2240/binning;
    angles          = 1:360;
    source_origin   = 109.83 + DistanceOffsetSample;
    origin_det      = DistanceSourceDetector - DistanceSourceOrigin;


    vol_geom = astra_create_vol_geom(dims(2),dims(1));
    proj_geom = astra_create_proj_geom(detector_geom,...
                                        det_width,...
                                        det_count,...
                                        angles,...
                                        DistanceSourceOrigin,...
                                        DistanceOriginDetector);

    % For CPU-based algorithms, a "projector" object specifies the projection
    % model used. In this case, we use the "strip" model.
    proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

    % Get the projection matrix as a Matlab sparse matrix.
    A = opTomo('line_fanflat', proj_geom,vol_geom);
    
    % Free memory
    astra_mex_projector('delete', proj_id);
 
end

