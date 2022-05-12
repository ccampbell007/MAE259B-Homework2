% MAE 259B Homework 2
% Chapter 6: Problem 3
% Integrated twist of a rod

% REFERENCES
% 
% [1] 	K. M. Jawed, "Conservative Force and Potential Energy," in Discrete Simulation of Slender Structures. 
% [2] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 5, Los Angeles: University of California Los Angeles, 2022. 
% [3] 	K. M. Jawed, parallel_transport, University of California Los Angeles, 2022.
% [4] 	K. M. Jawed, signed_angle, University of California Los Angeles, 2022.
% Notation Order - bottom exponent_topexponent_timestep
clc; clear all; close all;
%% Given

N = 5;              % Nodes
dt = .05;           % Seconds

% DOF Vector
q1 = [2e-2,0,0,0,6.306e-3,1.898e-2,0,0,-1.602e-2,1.197e-2,0,0,-1.641e-2,-1.143e-2,0,0,5.673e-3,-1.918e-2,0];                                                
q2 = [2e-2,0,0,0,6.306e-3,1.898e-2,0,-2.006e-1,-1.463e-2,1.281e-2,-8.462e-3,2.191e-1,-1.441e-2,-8.443e-3,-1.827e-2,4.726e-1,7.321e-3,-1.712e-2,-1.796e-2];
q3 = [2e-2,0,0,0,6.306e-3,1.898e-2,0,-3.908e-1,-1.119e-2,1.474e-2,-1.497e-2,3.572e-1,-9.813e-3,3.396e-4,-3.337e-2,9.616e-1,1.117e-2,-9.421e-3,-3.688e-2];

% Nodes
x1_1 = q1([1 2 3]);
x2_1 = q1([5 6 7]);
x3_1 = q1([9 10 11]);
x4_1 = q1([13 14 15]);
x5_1=  q1([17 18 19]);

x1_2 = q2([1 2 3]);
x2_2 = q2([5 6 7]);
x3_2 = q2([9 10 11]);
x4_2 = q2([13 14 15]);
x5_2 =  q2([17 18 19]);

x1_3 =  q3([1 2 3]);
x2_3 =  q3([5 6 7]);
x3_3 =  q3([9 10 11]);
x4_3 =  q3([13 14 15]);
x5_3 =  q3([17 18 19]);

% Edges
e1_1 = x2_1 - x1_1;
e2_1 = x3_1 - x2_1;
e3_1 = x4_1 - x3_1;
e4_1 = x5_1 - x4_1;

e1_2 = x2_2 - x1_2;
e2_2 = x3_2 - x2_2;
e3_2 = x4_2 - x3_2;
e4_2 = x5_2 - x4_2;

e1_3 = x2_3 - x1_3;
e2_3 = x3_3 - x2_3;
e3_3 = x4_3 - x3_3;
e4_3 = x5_3 - x4_3;

% Tangents
t1_1 = e1_1/norm(e1_1);
t2_1 = e2_1/norm(e2_1);
t3_1 = e3_1/norm(e3_1);
t4_1 = e4_1/norm(e4_1);

t1_2 = e1_2/norm(e1_2);
t2_2 = e2_2/norm(e2_2);
t3_2 = e3_2/norm(e3_2);
t4_2 = e4_2/norm(e4_2);

t1_3 = e1_3/norm(e1_3);
t2_3 = e2_3/norm(e2_3);
t3_3 = e3_3/norm(e3_3);
t4_3 = e4_3/norm(e4_3);

% Space-Parallel Reference Frame Director
u1 = [-8.11e-1,-5.851e-1,0];
u2 = parallel_transport(u1, t1_1, t2_1);
u3 = parallel_transport(u2, t1_1, t2_1);
u4 = parallel_transport(u3, t1_1, t2_1);

% Time Parallel Reference Frame
a1_1_1 = u1;
a1_2_1 = u2;
a1_3_1 = u3;
a1_4_1 = u4;

a1_1_2 = parallel_transport(a1_1_1, t1_1, t1_2);
a1_2_2 = parallel_transport(a1_2_1, t2_1, t2_2);
a1_3_2 = parallel_transport(a1_3_1, t3_1, t3_2);
a1_4_2 = parallel_transport(a1_4_1, t4_1, t4_2);

a1_1_3 = parallel_transport(a1_1_2, t1_2, t1_3);
a1_2_3 = parallel_transport(a1_2_2, t2_2, t2_3);
a1_3_3 = parallel_transport(a1_3_2, t3_2, t3_3);
a1_4_3 = parallel_transport(a1_4_2, t4_2, t4_3);

% Time parallel reference frame a2
v1_1 = cross(t1_1,u1);
v2_1 = cross(t2_1,u2);
v3_1 = cross(t3_1,u3);
v4_1 = cross(t4_1,u4);

a2_1_1 = v1_1;
a2_2_1 = v2_1;
a2_3_1 = v3_1;
a2_4_1 = v4_1;

v1_2 = cross(t1_2,u1);
v2_2 = cross(t2_2,u2);
v3_2 = cross(t3_2,u3);
v4_2 = cross(t4_2,u4);

a2_1_2 = v1_2;
a2_2_2 = v2_2;
a2_3_2 = v3_2;
a2_4_2 = v4_2;

v1_3 = cross(t1_3,u1);
v2_3 = cross(t2_3,u2);
v3_3 = cross(t3_3,u3);
v4_3 = cross(t4_3,u4);

a2_1_3 = v1_3;
a2_2_3 = v2_3;
a2_3_3 = v3_3;
a2_4_3 = v4_3;
% Twist Angle
th1_1 = q1([4]);
th2_1 = q1([8]);
th3_1 = q1([12]);
th4_1 = q1([16]);
th5_1=  q1([19]);

th1_2 = q2([4]);
th2_2 = q2([8]);
th3_2 = q2([12]);
th4_2 = q2([16]);
th5_2=  q2([19]);

th1_3 = q3([4]);
th2_3 = q3([8]);
th3_3 = q3([12]);
th4_3 = q3([16]);
th5_3 = q3([19]);

% Material Frame directors t=0
m1_1_1 = cos(th1_1)*a1_1_1 + sin(th1_1)*a2_1_1;
m2_1_1 = -sin(th1_1)*a1_1_1 + cos(th1_1)*a2_1_1;

m1_2_1 = cos(th2_1)*a1_2_1 + sin(th2_1)*a2_2_1;
m2_2_1 = -sin(th2_1)*a1_2_1 + cos(th2_1)*a2_2_1;

m1_3_1 = cos(th3_1)*a1_3_1 + sin(th3_1)*a2_3_1;
m2_3_1 = -sin(th3_1)*a1_3_1 + cos(th3_1)*a2_3_1;

m1_4_1 = cos(th4_1)*a1_4_1 + sin(th4_1)*a2_4_1;
m2_4_1 = -sin(th4_1)*a1_4_1 + cos(th4_1)*a2_4_1;

% Material Frame directors t=.05
m1_1_2 = cos(th1_2)*a1_1_2 + sin(th1_2)*a2_1_2;
m2_1_2 = -sin(th1_2)*a1_1_2 + cos(th1_2)*a2_1_2;

m1_2_2 = cos(th2_2)*a1_2_2 + sin(th2_2)*a2_2_2;
m2_2_2 = -sin(th2_2)*a1_2_2 + cos(th2_2)*a2_2_2;

m1_3_2 = cos(th3_2)*a1_3_2 + sin(th3_2)*a2_3_2;
m2_3_2 = -sin(th3_2)*a1_3_2 + cos(th3_2)*a2_3_2;

m1_4_2 = cos(th4_2)*a1_4_2 + sin(th4_2)*a2_4_2;
m2_4_2 = -sin(th4_2)*a1_4_2 + cos(th4_2)*a2_4_2;

% Material Frame directors t=.1
m1_1_3 = cos(th1_3)*a1_1_3 + sin(th1_3)*a2_1_3;
m2_1_3 = -sin(th1_3)*a1_1_3 + cos(th1_3)*a2_1_3;

m1_2_3 = cos(th2_3)*a1_2_3 + sin(th2_3)*a2_2_3;
m2_2_3 = -sin(th2_3)*a1_2_3 + cos(th2_3)*a2_2_3;

m1_3_3 = cos(th3_3)*a1_3_3 + sin(th3_3)*a2_3_3;
m2_3_3 = -sin(th3_3)*a1_3_3 + cos(th3_3)*a2_3_3;

m1_4_3 = cos(th4_3)*a1_4_3 + sin(th4_3)*a2_4_3;
m2_4_3 = -sin(th4_3)*a1_4_3 + cos(th4_3)*a2_4_3;

% Reference Twist at t=0 (Time shifted from textbook)
dm2_2 = 0;
dm3_2 = 0;
dm4_2 = 0;

% Reference Twist at t=.05 (Time shifted from textbook)
dm2_3 = computeReferenceTwist(a1_1_2,a1_2_2,t1_2,t2_2);
dm3_3 = computeReferenceTwist(a1_2_2,a1_3_2,t2_2,t3_2);
dm4_3 = computeReferenceTwist(a1_3_2,a1_4_2,t3_2,t4_2);


% Reference Twist at t=.1
dm2_4 = computeReferenceTwist(a1_1_3,a1_2_3,t1_3,t2_3);
dm3_4 = computeReferenceTwist(a1_2_3,a1_3_3,t2_3,t3_3);
dm4_4 = computeReferenceTwist(a1_3_3,a1_4_3,t3_3,t4_3);

% Integrated Discrete Twist t=0
tk_2_2 = 0;
tk_3_2 = 0;
tk_4_2 = 0;

% Integrated Discrete Twist t=.05
tk_2_3 = (th2_2 - th1_2) + dm2_3;
tk_3_3 = (th3_2 - th2_2) + dm3_3;
tk_4_3 = (th4_2 - th3_2) + dm4_3;

% Integrated Discrete Twist t=.1
tk_2_4 = (th2_3 - th1_3) + dm2_4;
tk_3_4 = (th3_3 - th2_3) + dm3_4;
tk_4_4 = (th4_3 - th3_3) + dm4_4;