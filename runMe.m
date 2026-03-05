% run-me file for Sheffield University ptychography reconstruction code.
% Ensure the data file is stored in the working directory, or load the data
% seperately and comment out the 'load' statement.
% Parameters are set for running WASP.
% For RAAR: try recon.beta = 0.85
% For rPIE: try recon.alpha = 0.1, recon.beta = 1
% For ePIE: usually recon.alpha = recon.beta = 1
% DM and ER do not require tuning parameters.
% HIVE (parallel WASP) requires the additional fields "recon.numWorkers"
% and "recon.subIters"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Citation for this data and code:                                        %
% Andrew. M. Maiden, Wenjie Mei and Peng Li,                              %
% "WASP: Weighted Average of Sequential Projections for ptychographic     %
% phase retrieval,"                                                       %
% XXX, pp. XX-XX (2024).                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath('algorithms');
addpath('tools');
addpath('ptychograms');

% change the filename here to load different datasests. 
data_filepath = 'ptychograms/OpticalPtychoData_HARDWARE_27-02-2026-16-29-24.mat';
load(data_filepath);

% Algorithm
algorithm_name = 'ePIE'; % 'WASP', 'RAAR', 'rPIE', 'ePIE', 'DM', 'ER'

% Option to create initial probe
initProbe = circleGenerator(600,3);

% dimention scaling and flipping
%expt.dps = flipdim(expt.dps,1); % diffraction patterns
%expt.dps = flipdim(expt.dps,2);
%expt.positions.x = 1.16*expt.positions.x; % Coordinates
%expt.positions.y = 1.06*expt.positions.y;

% set the reconstruction parameters
recon.iters      = 500;
recon.alpha      = 2;          
recon.beta       = 1;
recon.gpu        = 1;
recon.upLimit    = 2;        

% run the algorithm
algorithm = str2func(algorithm_name);
[obj, probe] = algorithm(expt, recon, initProbe);

% Save intensity and phase as two colour-mapped images without using JVM/figures
obj_intensity = abs(obj);
obj_phase = angle(obj);
probe_intensity = abs(probe);
probe_phase = angle(probe);

nmap = 256;
cmap = gray(nmap);
obj_rgb_in = mapToRGB(obj_intensity, cmap, nmap);
probe_rgb_in = mapToRGB(probe_intensity, cmap, nmap);
cmap = hsv(nmap);
obj_rgb_ph = mapToRGB(obj_phase, cmap, nmap);
probe_rgb_ph = mapToRGB(probe_phase, cmap, nmap);


% Add 1mm Scale Bar AND Label (Headless / No-JVM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-calculate pixel pitch (dx) in the sample plane
M_dps = size(expt.dps, 1);
N_dps = size(expt.dps, 2);
dx_sample = expt.wavelength * expt.cameraLength ./ ([M_dps, N_dps] * expt.cameraPixelPitch);
pixel_size_m = mean(dx_sample); 

% Configuration
scale_bar_m = 1e-3;         % 1 mm
bar_height_px    = 6;       % Thickness of the bar

% Calculate width in pixels
bar_width_px = round(scale_bar_m / pixel_size_m);

% Text Mask Generation ("1 mm")
char_1 = [0 1 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0; 0 1 1 1 0];
char_m = [0 0 0 0 0; 0 0 0 0 0; 1 1 1 1 0; 1 0 1 0 1; 1 0 1 0 1; 1 0 1 0 1; 1 0 1 0 1];
char_space = zeros(7, 2);

% Combine into one text block: "1" + space + "m" + space + "m"
text_mask = [char_1, char_space, char_m, char_space, char_m];

% Resize text mask to be larger (scale up by factor of 2 or 3 for visibility)
scale_text = 2; 
text_mask = kron(text_mask, ones(scale_text)); % Nearest-neighbor upscaling

[t_h, t_w] = size(text_mask);

% --- 4. Burn into Images ---
% Put images in a cell array to process them independently
images = {obj_rgb_in, obj_rgb_ph, probe_rgb_in, probe_rgb_ph};

for i = 1:length(images)
    img = images{i};
    [img_h, img_w, ~] = size(img); % Get dimensions for THIS specific image
    padding = 20; 

    % Coordinates for the Bar
    bar_r_start = img_h - padding - bar_height_px;
    bar_r_end   = img_h - padding;
    bar_c_start = img_w - padding - bar_width_px;
    bar_c_end   = img_w - padding;

    % Coordinates for the Text (Centered above the bar)
    text_r_start = bar_r_start - t_h - 5; % 5 pixels gap above bar
    text_r_end   = text_r_start + t_h - 1;
    
    % Center text horizontally relative to bar
    center_offset = floor((bar_width_px - t_w) / 2);
    text_c_start = bar_c_start + center_offset;
    text_c_end   = text_c_start + t_w - 1;
    
    % Only draw if the image is actually large enough to fit the bar and text
    if bar_r_start > 0 && bar_c_start > 0 && text_r_start > 0 && text_c_start > 0
        
        % --- ADAPTIVE COLOR LOGIC ---
        % Sample the image area where the scale bar will be placed
        bg_region = img(bar_r_start:bar_r_end, bar_c_start:bar_c_end, :);
        avg_brightness = mean(bg_region(:));
        
        % Determine image data type to set correct max/threshold values
        if isa(img, 'uint8')
            white_val = 255;
            black_val = 0;
            threshold = 127; % Midpoint of 0-255
        else
            % Assuming double or single precision (0.0 to 1.0)
            white_val = 1.0;
            black_val = 0.0;
            threshold = 0.5; % Midpoint of 0.0-1.0
        end
        
        % Decide the color based on the background brightness
        if avg_brightness < threshold
            draw_color = white_val; % Background is dark -> use white
        else
            draw_color = black_val; % Background is light -> use black
        end
        % ----------------------------

        for c = 1:3 % Loop R, G, B channels
            layer = img(:,:,c);
            
            % 1. Draw Bar
            layer(bar_r_start:bar_r_end, bar_c_start:bar_c_end) = draw_color;
            
            % 2. Draw Text (Only where mask is 1)
            roi = layer(text_r_start:text_r_end, text_c_start:text_c_end);
            roi(text_mask == 1) = draw_color;
            layer(text_r_start:text_r_end, text_c_start:text_c_end) = roi;
            
            img(:,:,c) = layer; % Put channel back
        end
    end
    images{i} = img; % Update cell array
end

% Unpack the modified images back to their original variables
obj_rgb_in   = images{1};
obj_rgb_ph   = images{2};
probe_rgb_in = images{3};
probe_rgb_ph = images{4};

fprintf('Added 1mm scale bar (%d px wide) and label to valid images.\n', bar_width_px);

if ~exist('results', 'dir')
    mkdir('results');
end

% Generate timestamp
now_str = datetime('now', 'Format', 'yyyy_MM_dd_HH_mm');
exp_id = regexp(data_filepath, '(?<=Data).*(?=\.mat)', 'match', 'once');

filename_base = sprintf('%s_%s_%d_%.2f_%.2f', now_str, algorithm_name, recon.iters, recon.alpha, recon.beta);

imwrite(obj_rgb_in, fullfile('results', sprintf('%s_obj_intensity%s.png', filename_base, exp_id)));
imwrite(obj_rgb_ph, fullfile('results', sprintf('%s_obj_phase%s.png', filename_base, exp_id)));
imwrite(probe_rgb_in, fullfile('results', sprintf('%s_probe_intensity%s.png', filename_base, exp_id)));
imwrite(probe_rgb_ph, fullfile('results', sprintf('%s_probe_phase%s.png', filename_base, exp_id)));