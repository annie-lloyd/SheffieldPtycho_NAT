    % helper to map a 2D array to an MxNx3 uint8 RGB image using the colormap
    function rgb = map_to_rgb(A, cmap, nmap)
        A(~isfinite(A)) = 0;                       % handle NaN/Inf
        Amin = min(A(:)); Amax = max(A(:));
        if Amin == Amax
            idx = ones(numel(A),1,'like',A);
        else
            idx = round((A(:)-Amin)./(Amax-Amin)*(nmap-1)) + 1;
            idx = min(max(idx,1),nmap);
        end
        rgb = reshape(uint8(255 * cmap(idx,:)), [size(A,1), size(A,2), 3]);
    end