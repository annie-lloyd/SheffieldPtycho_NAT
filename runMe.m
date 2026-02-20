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
load('OpticalPtychoDataExp11-02-2026-14-02-32.mat');

% Algorithm
algorithm_name = 'WASP'; % 'WASP', 'RAAR', 'rPIE', 'ePIE', 'DM', 'ER'

% Option to create initial probe
%initProbe = circleGenerator(512,60);

% dimention scaling and flipping
%expt.dps = flipdim(expt.dps,1); % diffraction patterns
%expt.dps = flipdim(expt.dps,2);
%expt.positions.x = 1.06*expt.positions.x; % Coordinates
%expt.positions.y = 1.06*expt.positions.y; 

% set the reconstruction parameters
recon.iters      = 5;
recon.gpu        = 1;            
recon.alpha      = 2;         
recon.beta       = 2;   
recon.upLimit    = Inf;
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

if ~exist('results', 'dir')
    mkdir('results');
end
imwrite(obj_rgb_in, fullfile('results', sprintf('%s_result_iter%04d_intensity.png', algorithm_name, recon.iters)));
imwrite(obj_rgb_ph, fullfile('results', sprintf('%s_result_iter%04d_phase.png', algorithm_name, recon.iters)));
imwrite(probe_rgb_in, fullfile('results', sprintf('%s_probe_iter%04d_intensity.png', algorithm_name, recon.iters)));
imwrite(probe_rgb_ph, fullfile('results', sprintf('%s_probe_iter%04d_phase.png', algorithm_name, recon.iters)));