%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Demo of CLIPPER+ max-clique point cloud registration 
% Thid demo is tested in Ubuntu 22.04, Matlab R2023a
%
% Download and uncompress the pointcloud datasets, and place them in the "data" folder
% Download link: https://drive.google.com/drive/folders/1-SshbPvfBeVXw3r7OazfwO0A1kwNSy23?usp=sharing
% Run this demo AFTER compiling CLIPPER+ with MATLAB binders.
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

% putative and inlier associations
fileID = fopen(strcat(adr_scan_pair, 'assoc_putative.txt'),'r');
assoc_putative = fscanf(fileID, '%d %d', [2, Inf]).';
fileID = fopen(strcat(adr_scan_pair, 'assoc_inliers.txt'),'r');
assoc_inliers = fscanf(fileID, '%d %d', [2, Inf]).';
fprintf('number of putative associations = %d\n', size(assoc_putative,1));
fprintf('number of inlier associations = %d\n', size(assoc_inliers,1));

% feature point coordinates in the pointclouds
fileID = fopen(strcat(adr_scan_pair, 'feat1.txt'),'r');
feat1 = fscanf(fileID, '%f %f %f', [3, Inf]).';
fileID = fopen(strcat(adr_scan_pair, 'feat2.txt'),'r');
feat2 = fscanf(fileID, '%f %f %f', [3, Inf]).';

% ground truth relative pose
fileID = fopen(strcat(adr_scan_pair, 'transf_gt.txt'),'r');
transf_gt = fscanf(fileID, '%f', [4, Inf]).';

% graph adjacency matrix
fileID = fopen(strcat(adr_scan_pair, 'adj.txt'),'r');
adj = fscanf(fileID, '%d', Inf);
siz_adj = sqrt(length(adj));
adj = reshape(adj, [siz_adj, siz_adj]);
fprintf('graph adjacency matrix size = %d\n', siz_adj);



























%% run Clipper+

fprintf('running clipper+ ...\n');
[omega_hat, clique_idx, certified, runtime] = clipperplus_clique_mex(adj);

fprintf('clipper+ done\n\n');


%% Display parameters

disp_figs = false; % general display flag, if 'false' no figs will show 
color_pt1 = [0, 0.4470, 0.7410]; % color of point cloud 1 
color_pt2 = [0.5412, 0.1686, 0.8863]; % color of point cloud 2
color_background = [1, 1, 1]; % backgroud color of point cloud figure

bench_file_format = '.ply';
if strcmp(bench_name, '7-Scenes') || strcmp(bench_name, 'Sun3D')
    bench_file_name = 'cloud_bin_';
end
if strcmp(bench_name, 'ETH')
    bench_file_name = 'Hokuyo_';
end

%% disply pointclouds

adr_pt = strcat(adr_scans, bench_file_name); % address to point cloud files 
fprintf('reading pointclouds...\n')
pt1_orig = pcread( strcat(adr_pt, num2str(idx1), bench_file_format) ); % first point cloud
pt2_orig = pcread( strcat(adr_pt, num2str(idx2), bench_file_format) ); % second point cloud

% add color
pt1_orig.Color = repmat(color_pt1,pt1_orig.Count,1);
pt2_orig.Color = repmat(color_pt2,pt2_orig.Count,1);

% display the pt clouds
figure;        
pcshow(pt1_orig, 'MarkerSize',3, 'BackgroundColor', color_background);
title('1st pointcloud');
view([0.2822 -66.9150])

figure;        
pcshow(pt2_orig, 'MarkerSize',3, 'BackgroundColor', color_background);
title('2nd pointcloud');
view([0.2822 -66.9150])

    

%% display aligned point clouds based on CLIPPER solution

% point cloud alignment
assoc_maxclq = assoc_putative(clique_idx,:); % max clique associations
feat1_maxclq = feat1(assoc_maxclq(:,1),:); % selected features
feat2_maxclq = feat2(assoc_maxclq(:,2),:); % selected features
[rot, trans] = arun(feat1_maxclq.', feat2_maxclq.'); % alignment rotation/translation from Arun's method

% check if correct solution found
rot_gt = transf_gt(1:3,1:3); % ground truth rotation
trans_gt = transf_gt(1:3,4); % ground truth translation
rot_err = norm( rotmat2vec3d(rot*rot_gt.') ); % rotation error in radians
trans_err = norm( trans - trans_gt); % translation error in meters

thresh_rot_err = 0.0873; % (in radian) 0.0873 radians is 5 degrees
thresh_trans_err = 2*0.05; % (in meters) 
if (rot_err < thresh_rot_err) && (trans_err < thresh_trans_err)
    fprintf('point clouds correctly registered by max clique solution\n')
else
    fprintf('point clouds wrongly registered by max clique solution\n')
end
fprintf('rotation error=%g radians;  translation error=%g meters\n', rot_err,trans_err);

% transform point cloud
transf = [rot, trans; [0 0 0 1]]; 
tform = rigidtform3d(transf); % matlab rigid transform
pt2_transf = pctransform(pt2_orig,tform); % transformed point cloud
feat2_maxclq_transf = transf * [feat2_maxclq.'; ones(1, size(feat2_maxclq,1))]; % transform feature points
feat2_maxclq_transf = feat2_maxclq_transf(1:3,:).';
residue_vec = sqrt(sum((feat1_maxclq - feat2_maxclq_transf).^2,2)); % distance between aligned features
fprintf('max distance error between aligned features = %g \n',max(residue_vec))

% display aligned pointclouds and aligned features (based on max clique solution)  
figure;
pt1_feat = pointCloud(feat1_maxclq, 'Color',[1,0,0]);
pt2_feat = pointCloud(feat2_maxclq_transf, 'Color',[0,1,0]);

hold on 
pcshow(pt1_orig, 'MarkerSize',30, 'BackgroundColor', color_background);
pcshow(pt2_transf, 'MarkerSize',30, 'BackgroundColor', color_background);
pcshow(pt1_feat, 'MarkerSize',500, 'BackgroundColor', color_background);
pcshow(pt2_feat, 'MarkerSize',500, 'BackgroundColor', color_background);
hold off
title('Aligned point clouds using CLIPPER ')        


