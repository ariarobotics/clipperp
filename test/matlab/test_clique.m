%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLIPPER+ clique test: finds a maximal clique in a given graph, and 
% certify if it's the maximum clique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

addpath(genpath('../../build/bindings/matlab')) % clipper+ mex file
%%

% graph adjacency matrix
adj = [ 0, 0, 1, 1, 1;
        0, 0, 1, 1, 1;
        1, 1, 0, 1, 1;
        1, 1, 1, 0, 1;
        1, 1, 1, 1, 0];   

% run clipper+ clique finding algorithm
[clique_size, clique, certificate, runtime] = clipperplus_clique_mex(adj)
