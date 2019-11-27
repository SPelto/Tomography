function A = createSysMat(binning,angles)
    
   %%%% Creates forward projector operator. %%%%
   % binning    = binningfactor
   % angles     = angles of the x-ray source
    
    % Setup of the x-ray device
    DistanceOffsetSample    = 275; % Bone
    DistanceSourceDetector  = 552.18;
    DistanceSourceOrigin    = 109.83 + DistanceOffsetSample;
    DistanceOriginDetector  = DistanceSourceDetector - DistanceSourceOrigin;
    pixelSize               = 0.050 * binning;
    M                       = DistanceSourceDetector / DistanceSourceOrigin;
    effPixelSize            = pixelSize / M;
    
    % Distance from source to origin specified in terms of effective pixel size
    DSO             = DistanceSourceOrigin / effPixelSize;
    % Distance from origin to detector specified in terms of effective pixel size
    DOD             = DistanceOriginDetector /effPixelSize;

    % Some necessary parameters for ASTRA
    detector_geom   = 'fanflat';
    det_count       = 2240/binning;
    angles          = deg2rad(angles);
    reconSize       = 2048/binning;


    % Create the geometries necessary for ASTRA
    vol_geom = astra_create_vol_geom(reconSize,reconSize);
    proj_geom = astra_create_proj_geom(detector_geom,...
                                        M,...
                                        det_count,...
                                        angles,...
                                        DSO,...
                                        DOD);

    % Get the projection matrix as a Matlab sparse matrix.
    A = opTomo('strip_fanflat', proj_geom,vol_geom);
    
    % Memory cleanup of mex files
    astra_mex_data2d('delete', vol_geom);
    astra_mex_data2d('delete', proj_geom);
end

