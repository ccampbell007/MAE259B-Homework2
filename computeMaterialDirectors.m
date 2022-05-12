function [m1,m2] = computeMaterialDirectors(d1,d2,theta);

% Calcualte Material Director Frame

ne = numel(theta);  % Number of Edges

% Initialize vectors
m1 = zeros(ne,3);   
m2 = zeros(ne,3);

% Calculate Material Directors
for c = 1:ne
    cs = cos(theta(c));
    ss = sin(theta(c));
    d1_l = d1(c,:);
    d2_l = d2(c,:);
    m1_l = cs * d1_l + ss * d2_l;
    m1(c,:) = m1_l / norm(m1_l);
    m2_l = - ss * d1_l + cs * d2_l;
    m2(c,:) = m2_l / norm(m2_l);
end