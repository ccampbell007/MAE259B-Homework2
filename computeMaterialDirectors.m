function [m1,m2] = computeMaterialDirectors(d1,d2,theta);
% Reference:
% [7] K. M. Jawed, computeMaterialDirectors, University of California Los Angeles, 2022. 


% Calcualte Material Director Frame

ne = numel(theta);  % Number of Edges

% Initialize vectors
m1 = zeros(ne,3);   
m2 = zeros(ne,3);

% Calculate Material Directors
for k = 1:ne
    d1_l = d1(k,:);
    d2_l = d2(k,:);
    m1_l = cos(theta(k)) * d1_l + sin(theta(k)) * d2_l;
    m1(k,:) = m1_l / norm(m1_l);
    m2_l = - sin(theta(k)) * d1_l + cos(theta(k)) * d2_l;
    m2(k,:) = m2_l / norm(m2_l);
end