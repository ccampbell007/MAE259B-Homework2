% MAE 259B Homework 2
% Chapter 7
% Integrated twist of a rod

% REFERENCES
% 
% [1] 	K. M. Jawed, "Conservative Force and Potential Energy," in Discrete Simulation of Slender Structures. 
% [2] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 12, Los Angeles: University of California Los Angeles, 2022. 
% [3] 	K. M. Jawed, MAE 259B - Spring 2022 - Lecture 13, Los Angeles: University of California Los Angeles, 2022.
% [4] 	K. M. Jawed, computekappa, University of California Los Angeles, 2022. 
% [5] 	K. M. Jawed, parallel_transport, University of California Los Angeles, 2022.
% [6] 	K. M. Jawed, signed_angle, University of California Los Angeles, 2022.
% [7] 	K. M. Jawed, computeMaterialDirectors, University of California Los Angeles, 2022. 
% [8] 	K. M. Jawed, computeReferenceTwist, University of California Los Angeles, 2022. 
% [9] 	K. M. Jawed, computeTangent, University of California Los Angeles, 2022. 
% [10] 	K. M. Jawed, computeTimeParallel, University of California Los Angeles, 2022. 
% [11] 	K. M. Jawed, crossMat, University of California Los Angeles, 2022.
% [12] 	K. M. Jawed, getkappa, University of California Los Angeles, 2022.
% [13] 	K. M. Jawed, getRefTwist, University of California Los Angeles, 2022.
% [14] 	K. M. Jawed, gradEb_hessEb, University of California Los Angeles, 2022.
% [15] 	K. M. Jawed, gradEs_hessEs, University of California Los Angeles, 2022.
% [16] 	K. M. Jawed, gradEt_hessEt, University of California Los Angeles, 2022.


clc; clear all; close all;
%% Establish Initial Conditions Prior To Using Newton's Law
% Need Edges, reference directors, material directors, tangents
% Parallel transport between nodes
%% Given Properties

% Rod
N = 50;                         % Number of nodes
l = .20;                        % Length (m)
Rn = .02;                       % Radius (m)
dtheta = 1/Rn*1/(N-1);          % Theta increment (rad)
dt = .01;                       % Time Increment (s)
RunTime = 5;                    % Run Time (s)
Nsteps = RunTime/dt;            % Number of steps
rho = 1000;                     % Density (kg/m^3)
r0 = .001;                      % Cross-sectional radius (m)
nu = 0.5;                       % Poission Ratio
Y = 10e6;                       % Elastic Modulus (Pa)
G = Y/(2*(1+nu));               % Shear Modulus (Pa)
g = [0,0,-9.81];                % Gravity (m/s^2)
I = pi/4*(r0^4);                % Mass Moment of Intertia (m^4)
J = pi/2*r0^4;                  % Second Polar Moment of Area
A = pi*(r0^2);                  % Beam Cross Sectional Area (m^2)

% Simplified Material Properties
EA = Y*A;                       
EI = Y*I;
GJ = G*J;

% Counting and mass variables
ndof = 4*N - 1;                            % Degrees of freedom
ne = N - 1;                                % Number of edges
dm = (pi*r0^2*l)*rho/ne;                   % Node mass
ScaleSolver = EI/l^2;                      % EI and rodlength relation from class

%% Fixed/Free nodes
% Fixed and free nodes
fixed_index = 1:7;
free_index = 8:ndof;

%% Node matrix
nodes = zeros(N,3);

% Incremental Angle
dTheta = (l/ne)/Rn;                        

% Initial Nodes (Given)
for c = 1:N
    nodes(c,1) = Rn*cos((c-1)*dTheta);
    nodes(c,2) = Rn*sin((c-1)*dTheta);
end

%% Initial Locations
% Initialize initial locations
x0 = zeros(ndof,1);
for c = 1:N
    x0(4*(c-1) + 1) = nodes(c,1); % x-coord
    x0(4*(c-1) + 2) = nodes(c,2); % y-cood
    x0(4*(c-1) + 3) = nodes(c,3); % z-coord
end

% Theta (4n element)
x0(4:4:end) = 0;
x = x0;           % New DOF

%% Initial Velocity
qdot = (x-x0)/dt; 

%% Define DOF Vector
% DOF vector (x,y,z,theta)
for k=1:N    
q(4*k-3,1) = x(k,1);         % X Coordinate
q(4*k-2,1) = x(2*k,1);         % Y Coordinate
q(4*k-1,1) = x(3*k,1);         % Z Coordinate
q(4*k,1) = 0;                % Angle of twist
end                 
q(end) = [];                % Remove last term

%% Mass Matrix
% Initialize matrix
m = zeros(ndof,1);

% Node Mass Values (rods)
for i = 1:N
    if i == 1 || i == N                 % First and last node
        m(4*(i-1)+1:4*(i-1)+3) = dm/2;
    else
        m(4*(i-1)+1:4*(i-1)+3) = dm;    % MIddle nodes
    end    
end

% Node Masses
for c = 1:ne
    m(4*c) = dm/2*r0^2;
end

%% Gravity
% Initialize Matrix
w = zeros(ndof,1);

% Assign gravity to x,y,z nodes
for i = 1:N
    w(4*(i-1)+1:4*(i-1)+3) = g;
end

Fg = m.*w;   % Weight

%% Reference Length
% Initialize Matrix
dist = zeros(ne,1);

% Distance between nodes
for c = 1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    dist(c) = norm(dx);
end

%% Veronoi Length
% Initialize Matrix
veronoiRefLen = zeros(N,1);

% Average edge lengths
for c = 1:N
    if c==1
        veronoiRefLen(c) = 0.5*dist(c);
    elseif c==N
        veronoiRefLen(c) = 0.5*dist(c-1);
    else
        veronoiRefLen(c) = 0.5*(dist(c-1) + dist(c));
    end
end

%% Reference Directors
% Initialize directors and tangent
d1 = zeros(ne,3);
d2 = zeros(ne,3);
tangent = zeros(ne,3);

% Differential distance and tangent
for c=1:ne
    dx = nodes(c+1,:) - nodes(c,:);
    tangent(c,:) = dx/norm(dx);
end

% Initial directors and tangents
t0 = tangent(1,:);
t1 = [0,0,-1];
d11 = cross(t0,t1);

% Space-Parallel Reference Frame Director
%u(1,:) = [0,0,1];

% For small d1 values, establish initial conditions
if abs(d11) < 1e-6
    t1 = [0 1 0];
    d11 = cross(t0,t1);
end

% Directors
d1(1,:) = d11 / norm(d11);
d21 = cross(t0, d1(1,:));
d2(1,:) = d21 / norm(d21);

% Directors for each node via parallel transport
for c = 2:ne
    t0 = tangent(c-1,:);
    t1 = tangent(c,:);
    d1_old = d1(c-1,:);
    d1_new = parallel_transport(d1_old, t0, t1);
    d1_new = d1_new / norm(d1_new);
    d1(c,:) = d1_new;
    d2_new = cross(t1, d1_new);
    d2(c,:) = d2_new / norm(d2_new);
end

% % Define Initial Edges
% for i = 1:N-1
%     e(i,:) = x(i+1,:) - x(i,:);
% 
%     % Tangents
% %    t(i,:) = e(i,:)/norm(e(i,:));
% end

%% Material Director
% for i = 1:N-1
% % Material Frame directors t=0
% %% Reference Frames
% % Time-Parallel and Space-Parallel Reference Frames
% a1(i,:) = u(i,:);
% u(i+1,:) = parallel_transport(u(i,:), tangent(1,:), tangent(2,:));
% 
% % Time parallel reference frame a2
% v(i,:) = cross(tangent(i,:),u(i,:));
% 
% % Time Parallel Reference Frame
% a2(i,:) = v(i,:);
% 
% end

%% Material Director
% Initial Material Directior
theta = x(4:4:end);
[m1,m2] = computeMaterialDirectors(d1,d2,theta);

% Reference Twist
% Initialize
refTwist = zeros(N,1);
% Initial Reference Twist
refTwist = getRefTwist(d1,tangent,refTwist);

% for i = 1:N-2
% m1(i,:) = cos(theta(i,:))*a1(i,:) + sin(theta(i+1,:))*a2(i,:);
% m2(i,:) = -sin(theta(i,:))*a1(i,:) + cos(theta(i,:))*a2(i,:);
% 
% end

%% Natural Curvature
kappa = getkappa(x,m1,m2);

%% Timstep Loop
Nsteps = round(RunTime/dt); % Number of timesteps
z_end = zeros(Nsteps,1);    % Value to store
ctime = 0;                  % Current Time
d1_old = d1;                % Re-assign director
d2_old = d2;                % Re-assign director
refTwist_old = refTwist;    % Re-assign reference twist
ii = 0;                     % Counter

% Loop for each timestep
for timestep = 1:Nsteps
    ii = ii + 1;             % Increment
    fprintf('t-%f\n',ctime); % Print Timestep

    % Newton's Law Function
    [x,d1,d2,refTwist] = DISCRETE_ELASTIC_RODS(q,qdot,tangent,dt,EI,EA,GJ,N,free_index,ScaleSolver,x0,d1_old,m,refTwist,Fg,veronoiRefLen,kappa,refTwist_old,dist);
      
   % Plot 3d profile
   x_coord = x(1:4:end);
   y_coord = x(2:4:end);
   z_coord = x(3:4:end);
   figure(1);
   clf();
   plot3(x_coord, y_coord, z_coord)
   axis equal
   xlabel('x Location (m)')
   ylabel('y Location (m)')
   zlabel('z Location (m)')
   
   
   % Store Z Value
   z_end(timestep) = x(end);
  
   % Update velocity, position, directors, reference twist
   qdot = (x-x0)/dt; 
   x0 = x; 
   d1_old = d1;
   d2_0ld = d2;
   refTwist_old = refTwist;

   % Increment Time
    ctime = ctime + dt;
end

 % Plot last node z coordinate
 figure(2);
 tarray = (1:Nsteps) + dt;
 plot(tarray,z_end,'b-*')
 xlabel('Timestep')
 ylabel('Tip displacement (m)')
title('Tip Displacement vs Time')
