%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Demo of CLIPPER+ max-clique  
% Thid demo is tested in Ubuntu 22.04, Matlab R2023a
%
% Download and uncompress the pointcloud datasets, and place them in the "data" folder
% Download link: https://drive.google.com/drive/folders/1-SshbPvfBeVXw3r7OazfwO0A1kwNSy23?usp=sharing
% Run this demo AFTER compiling CLIPPER+ with MATLAB binders.
% 
%
% This Demo was tested in Ubuntu 22.04, and MATLAB R2023a
% (C) Kaveh Fathian, Tyler Summers, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("./helpers")

% add clipper+ library & mex file to path
addpath("../build/");
addpath("../build/bindings/matlab");

adr_data = './data'; % address to dataset 


%% set benchmark and sequence

% choose dataset 7 sequence
bench_name = '7-Scenes';
bench_sequence = 'kitchen';

% bench_name = 'Sun3D';
% bench_sequence = 'sun3d-home_at-home_at_scan1_2013_jan_1'; % {'sun3d-home_at-home_at_scan1_2013_jan_1', 'sun3d-hotel_uc-scan3', 'sun3d-mit_76_studyroom-76-1studyroom2'};

% bench_name = 'ETH';
% bench_sequence = 'gazebo_summer'; % {'gazebo_summer', 'wood_autmn'};


%% matched pointclouds

% index of matched pointclouds. change to see different matches (matched 
% pairs must exist in the dataset folder)
idx1 = 0; 
idx2 = 1;


%% dataset address

adr_scans = strcat(adr_data, '/', bench_name, '/', bench_sequence, '/');        
folder_name = strcat('scans_', num2str(idx1), '_', num2str(idx2));
adr_scan_pair = strcat(adr_scans, folder_name, '/');


%% import data

fprintf('loading scans %d and %d of in %s dataset.\n', idx1, idx2, bench_name);

% graph adjacency matrix
fileID = fopen(strcat(adr_scan_pair, 'adj.txt'),'r');
adj = fscanf(fileID, '%d', Inf);
siz_adj = sqrt(length(adj));
adj = reshape(adj, [siz_adj, siz_adj]);
fprintf('graph adjacency matrix size = %d\n', siz_adj);

% ground truth maximum clique and its index (in the adjacency matrix)
fileID = fopen(strcat(adr_scan_pair, 'omega_gt.txt'),'r');
omega_gt = fscanf(fileID, '%d', Inf).';
fileID = fopen(strcat(adr_scan_pair, 'omega_gt_idx.txt'),'r');
omega_gt_idx = fscanf(fileID, '%d', Inf);
fprintf('ground truth clique size = %d\n', omega_gt);


%% run Clipper+

fprintf('running clipper+ ...\n');
[omega_hat, clique_idx, certified, runtime] = clipperplus_clique_mex(adj);

fprintf('clipper+ done\n\n');

omega_ratio = omega_hat/omega_gt;
omega_ratio2 = numel(intersect(omega_gt_idx,clique_idx))/omega_gt;


fprintf('clipper+ found clique of size %d.\n', omega_hat);
fprintf('ground truth clique size is %d.\n', omega_gt);



