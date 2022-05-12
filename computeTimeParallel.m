function [d1,d2] = computeTimeParallel(d1_old,x0,x,N)
% Computes new directors as time progresses
nv = (numel(x)+1) / 4;
ne = N -1;

% Calculate tangents
tangent0 = zeros(ne,3);
for c = 1:ne
    dx = x0( 4*c+1:4*c+3 ) - x0( 4*(c-1)+1:4*(c-1)+3 );
    tangent0(c,:) = dx / norm(dx);
end

tangent = zeros(ne,3);
for c = 1:ne
    dx = x( 4*c+1:4*c+3 ) - x( 4*(c-1)+1:4*(c-1)+3 );
    tangent(c,:) = dx / norm(dx);
end

% Initialize Matrices
d1 = zeros(ne,3);
d2 = zeros(ne,3);

% Calculate new directors
for c = 1:ne
    d1_l = parallel_transport(d1_old(c,:), tangent0(c,:), tangent(c,:));
    d1_l = (d1_l - dot(d1_l,tangent(c,:)) * tangent(c,:)); 
    d1_l = d1_l / norm(d1_l);
    d1(c,:) = d1_l;
    d2_l = cross( tangent(c,:), d1_l);
    d2(c,:) = d2_l / norm(d2_l);
end
end