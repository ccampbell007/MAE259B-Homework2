function [x,d1,d2,refTwist] = DISCRETE_ELASTIC_RODS(q,qdot,tangent,dt,EI,EA,GJ,N,free_index,ScaleSolver,x0,d1_old,m,refTwist,Fg,veronoiRefLen,kappa,refTwist_old,dist)
% Run Netwon's Law Simulation fo Discrete Twist

q0 = q;                     % Initial Guess
x = x0;                     % Initial Guess
n = 0;                      % Counting Variable
tol = 1e-3;                 % Tolerance
err = tol*ScaleSolver*10;   % Set initial error value
maximum_iter = 100;         % Maximum iteration

% Newton's Law Loop
while err > tol*ScaleSolver
% Increment
n = n + 1;
% Update Reference Directors and tangent
[d1,d2] = computeTimeParallel(d1_old,x0,x,N);
tangent = computeTangent(x,N);
 
%     for i = 1:N-1
%     e(i,:) = q(i+1) - q(i);
%     tangent(n+1,:) = e(i,:)/norm(e(i,:));
%     th(i) = q(4*i,n);
%  
%     % Compute Reference Frame % Check this
%     a1(i,:) = parallel_transport(a1(i,:),tangent(n,:),tangent(n+1,:));
%     a2(i,:) = parallel_transport(a2(i,:),tangent(n,:),tangent(n+1,:)); % *Questionable*
%     
%     % Compute material frame (m1,m2,t)@t+dt
%     m1(i,:) = cos(th(:,i)).*a1(i,:) + sin(th(:,i)).*a2(i,:);
%     m2(i,:) = -sin(th(:,i)).*a1(i,:) + cos(th(:,1)).*a2(i,:);
%     
%     % Compute Referenec Twist * Double Check*
%     dm(i,:) = computeReferenceTwist(a1(i,:),a2(i,:),tangent(n,:),tangent(n+1,:));
%     end

% Reference Twist
refTwist = getRefTwist(d1,tangent,refTwist_old);

%Material Director
theta = x(4:4:end);
[m1,m2] = computeMaterialDirectors(d1,d2,theta);
    
% Get Forces
% Calculate Elastic Energy Terms
[Fs,Ft,Fb,Js,Jt,Jb] = GRAD_HESS_ELASTIC_ROD(x,m1,m2,refTwist,EI,EA,GJ,N,kappa,theta,veronoiRefLen,dist);

% Calculate F and J
F = m.*(x-x0)/dt^2 - m.*qdot/dt - (Fb+Ft+Fs+Fg);
J = diag(m)/dt^2 - (Jb+Jt+Js);

% Calculate location change
F_free = F(free_index);
J_free = J(free_index,free_index);
dx = J_free\F_free;

% Update location
x(free_index) = x(free_index) - dx;

% Check error
err = sum(abs(F_free));
fprintf('err - %f, iter - %d\n',err,n);
if n > maximum_iter
        fprintf('Error\n');
        return
    end
end
end

%  % Initial Position and Velocity
%     q(free_index) = x(free_index);    
%     q = q0;                     % DOF
%     qdot = (q - q0)/dt;            % Velocity
%     % Compute Reference Frame (a1,a2,t)@t+dt using q(t+dt)
%     
%     % Calculate tangent and angle at new timestep
%     for i = 1:N
%     e(i,:) = q(i+1) - q(i);
%     tangent(n+1,:) = e(i,:)/norm(e(i,:));
%     if i == N
%     else
%     th(i) = q(4*i);
%     end
%     end
%     
%     for i = 1:N-1
%     % Compute Reference Frame % Check this
%     a1(i,:) = parallel_transport(a1(i,:),tangent(n,:),tangent(n+1,:));
%     a2(i,:) = parallel_transport(a2(i,:),tangent(n,:),tangent(n+1,:)); % *Questionable*
%     end
    %% Extra Code

    % Define q(t+dt)    
    
%     % Compute Reference Frame (a1,a2,t)@t+dt using q(t+dt)
%     
%     % Calculate tangent and angle at new timestep
% 
%     for i = 1:N-1
%     % Compute Reference Frame % Check this
%     a1(i,:) = parallel_transport(a1(i,:),t(n,:),t(n+1,:));
%     a2(i,:) = parallel_transport(a2(i,:),t(n,:),t(n+1,:)); % *Questionable*
% 
%     % Compute material frame (m1,m2,t)@t+dt
%     m1(i,:) = cos(th(:,i)).*a1(i,:) + sin(th(:,i)).*a2(i,:);
%     m2(i,:) = -sin(th(:,i)).*a1(i,:) + cos(th(:,1)).*a2(i,:);
%     end
%     
%     % Compute Referenec Twist * Double Check*
%     dm(i,:) = computeReferenceTwist(a1(i,:),a2(i,:),t(n,:),t(n+1,:));
%     
%     % Interia Term
%     f = M/dt * ((q - q0)/dt - qdot);
%     J = M/dt^2;
%       
%     % Weight Term
%     f = f - W;
%     
%     % Elastic Energy Terms *Add global variables*
%     [f,J] = GRAD_HESS_ELASTIC_ROD(q,m1,m2,a1,a2,dm,EI,EA,GJ,N);
% 
%     % Update position
%     f_free = f(free_index);
%     J_free = J(free_index, free_index);
%     q_free = q_free - J_free \ f_free;                    
%     
%     err = sum(abs(f_free));         % Evaluate error to continue N-R Method
%     
%     q(free_index) = q_free;
% 
%     n=n+1;
%     end
%     
% end
