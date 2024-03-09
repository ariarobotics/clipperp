%% test CLIPPER+ maximal clique on DIMACS graph dataset
%
% THE SECOND DIMACS IMPLEMENTAION CHALENGE (1992-1993) 
% https://iridia.ulb.ac.be/~fmascia/maximum_clique/DIMACS-benchmark
% 

clear all;
close all;
clc;

addpath(genpath('../../build/bindings/matlab')) % clipper+ mex file

%% load dataset
 
Ms = load('DIMACS9.mat');

names = {'M_p_hat300_1'; 'M_p_hat300_2'; 'M_C1259'; 'M_C2509'; ...
    'M_keller4'; 'M_p_hat300_2'; 'M_brock200_4'; ...
    'M_brock200_2'; 'M_gen200_p09_44'; 'M_gen200_p09_55'}; % graphs in DIMACS dataset
omega_gt = [8; 25; 34; 44; 11; 25; 17; 12; 44; 55]; % ground truth clique size


%% 

num_graphs = length(names); % number of graphs

omega_ratio = zeros(num_graphs, 1);
runtime = zeros(num_graphs, 1);
density = zeros(num_graphs, 1);
certificates = zeros(num_graphs, 1);

for graph = 1 : num_graphs 
    M = Ms.(names{graph}); % matrix M for this benchmark
    n = length(M); % nodes
    adj = M - eye(n); % adjacency matrix
    
    fprintf('processing graph %s ...\n', names{graph});
    
    density(graph) = nnz(adj) / (n*(n-1));
    fprintf('graph density = %g \n', density(graph));
    
    fprintf('benchmarking clipper+_clique\n');
    [omega_hat, clique_idx, certificate, t] = clipperplus_clique_mex(adj);
    omega_ratio(graph, 1) = omega_hat/omega_gt(graph); % max clique ratio at this run
    certificates(graph,1) = certificate;
    runtime(graph, 1) = t; 
    
    fprintf('\n\n');
end % end for


%% display results

or_clipper_plus = omega_ratio(:,1);

th = table(names, density, omega_gt, ...
    or_clipper_plus, runtime, certificates)










