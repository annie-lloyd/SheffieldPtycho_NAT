function grid = circleGenerator(image_size,r)

% circleGenerator generates a binary grid with a filled circle of radius r
% centered in the grid.
% Inputs:
%   image_size - Size of the square grid (e.g., 512 for a 512x512 grid)
%   r          - Radius of the circle
% Output:
%   grid       - Binary grid with a filled circle

    % Define grid dimensions
    width = image_size;
    height = image_size;

    % Create a meshgrid of coordinates
    % 1:width corresponds to columns (x), 1:height to rows (y)
    [x, y] = meshgrid(1:width, 1:height);

    % Calculate the center of the image
    % For a 512x512 image, the exact center is at 256.5
    cx = (width + 1) / 2;
    cy = (height + 1) / 2;

    % Calculate the squared distance from the center for every point
    % Using (x-cx)^2 + (y-cy)^2
    dist_sq = (x - cx).^2 + (y - cy).^2;

    % Generate the binary grid
    % Points where distance <= radius are 1, others are 0
    grid = double(dist_sq <= r^2);

end