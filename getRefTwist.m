function refTwist = getRefTwist(d1,tangent, refTwist)
% Reference:
% [13] 	K. M. Jawed, getRefTwist, University of California Los Angeles, 2022.


% Compute Reference Twist

% Initialize
[ne, ~] = size(d1);
refTwist = zeros(ne+1,1);

for c =2:ne
    % Parameters for computeReferenceTwist
    u0 = d1(c-1,:);
    u1 = d1(c,:);
    t0 = tangent(c-1,:);
    t1 = tangent(c,:);
    
    % Reference Twist Function
    reftwist(c) = computeReferenceTwist(u0,u1,t0,t1);
end

end