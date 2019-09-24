close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
%%%v(k,A)=max u +beta*E[v(k',A')]


%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
Pi_h=0.977;
Pi_l=0.926;
Pi = [Pi_h 1-Pi_h;1-Pi_l Pi_l]
A_h = 1.1;
A_l = 0.678;
A = [A_h A_l];

cons_h =  A(1)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat';
cons_l =  A(2)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % return function
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat_h = ret_h + beta * repmat(Pi(1,:)*v_guess, [num_k 1]);
    value_mat_l = ret_l + beta * repmat(Pi(2,:)*v_guess, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn = [vfn_h'; vfn_l'];
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

g1 = k(pol_indx_h); % policy function
g2 = k(pol_indx_l); % policy function

plot(k,vfn)
figure
plot(k,g1,g2)
plot(k,g2)

%%%simulatoin


