% MAE 259B Homework 2
% Chapter 6 Problem 1
% Integrated twist of a rod

% REFERENCES
% 
% [1] 	K. M. Jawed, "Conservative Force and Potential Energy," in Discrete Simulation of Slender Structures. 
% [2] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 5, Los Angeles: University of California Los Angeles, 2022. 
% [3] 	K. M. Jawed, parallel_transport, University of California Los Angeles, 2022.
% [4] 	K. M. Jawed, signed_angle, University of California Los Angeles, 2022.

clc; clear all; close all;
%% Given

% Nodes
x1 = [0.00, 0.00, 0.00];
x2 = [0.50, 0.00, 0.00];
x3 = [0.75, 0.25, 0.00];
x4 = [0.75, 0.50, 0.25];

% Material Directors
m1 = [0.00, 0.00, 1.00];
m2 = [0.00, 0.00, 1.00];
m3 = [0.00, -1/(2)^(0.5),1/(2)^(0.5)];

% Compute Vectors
u1 = [0,0,1];
v1 = [0,1,0];
t1 = [1,0,0];

u2 = [0,0,1];
v2 = [0,1,0];

% Edges
e1 = x2 - x1;
e2 = x3 - x2;
e3 = x4 - x3;

% Tangents
t2 = e2/norm(e2);
t3 = e3/norm(e3);

% Compute parallel transport
u2 = parallel_transport(u1, t1, t2);
u3 = parallel_transport(u2, t2, t3);
P2 = parallel_transport(m1, t1, t2);
P3 = parallel_transport(m2, t2, t3);

% Compute discrete angle of twist
tau2 = signedAngle(P2,m2,t2)
tau3 = signedAngle(P3,m3,t3)