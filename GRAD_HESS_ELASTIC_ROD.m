function [Fs,Ft,Fb,Js,Jt,Jb] = GRAD_HESS_ELASTIC_ROD(x,m1,m2,refTwist,EI,EA,GJ,N,kappa,theta,veronoiRefLen,dist)
% Gradient and Hessian of Elastic Energy in a Rod
% Initialize Vectors
ndof = numel(x);
Fs = zeros(ndof, 1);
Js = zeros(ndof,ndof);
Ft = zeros(ndof, 1);
Jt = zeros(ndof,ndof);
Fb = zeros(ndof, 1);
Jb = zeros(ndof,ndof);

% Stretching Elastic Energy
for k = 1:N-1
    node1 = [x(4*k-3), x(4*k-2),x(4*k-1)];
    node2 = [x(4*k+1), x(4*k+2),x(4*k+3)];
    ind = [4*k-3,4*k-2,4*k-1,4*k+1,4*k+2,4*k+3]; 
    [dF,dJ] = gradEs_hessEs(node1,node2,dist(k),EA);
    Fs(ind) = Fs(ind) - dF;
    Js(ind,ind) = Js(ind,ind) - dJ;
end

% Torsion Elastic Energy
for k=2:N-1
    node1 = [x(4*k-7), x(4*k-6),x(4*k-5)];
    node2 = [x(4*k-3), x(4*k-2),x(4*k-1)];
    node3 = [x(4*k+1), x(4*k+2),x(4*k+3)];
    theta_e = x(4*k-4);
    theta_f = x(4*k);
    ind = [4*k-7,4*k-6,4*k-5,4*k-4,4*k-3,4*k-2,4*k-1,4*k,4*k+1,4*k+2,4*k+3];
    [dF,dJ] = gradEt_hessEt(node1,node2,node3,theta_e,theta_f,refTwist(k),veronoiRefLen(k),GJ);
    Ft(ind) = Ft(ind) - dF;
    Jt(ind,ind) = Jt(ind,ind) - dJ;
end

for k=2:N-1
    node1 = [x(4*k-7), x(4*k-6),x(4*k-5)];
    node2 = [x(4*k-3), x(4*k-2),x(4*k-1)];
    node3 = [x(4*k+1), x(4*k+2),x(4*k+3)];
    m1e = m1(k-1,:);
    m2e = m2(k-1,:);
    m1f = m1(k,:);
    m2f = m2(k,:);
    ind = [4*k-7,4*k-6,4*k-5,4*k-4,4*k-3,4*k-2,4*k-1,4*k,4*k+1,4*k+2,4*k+3];
    [dF,dJ] = gradEb_hessEb(node1,node2,node3,m1e,m2e,m1f,m2f,kappa(k,:),veronoiRefLen(k),EI);
    Fb(ind) = Fb(ind) - dF;
    Jb(ind,ind) = Jb(ind,ind) - dJ;
end
end

%for k = 1:N
% x(k,1) = q(4*k-3,1);
% x(k,2) = q(4*k-2,1);
% x(k,3) = q(4*k-1,1);

% % Edge length
% for k = 0:N-2
%     e(k+1,:) = -x(4*k+1:4*k+3) + x(4*k + 5:4*k+7);
% end
% 
% % Veronoi Length
% for k = 2:N-1
%     l_k(k,:) = (norm(e(k-1,:))+norm(e(k,:)))/2;
% end