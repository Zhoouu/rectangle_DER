function [gradE, hessE, potE] = computeContactPotentialP2P(x1s, x2s, K1, h2, state, contact_stiffness, scale)
% Calculates P2P potential gradient, Hessian, and Energy. (Analytic)

    gradE_local = zeros(6, 1); % Local gradient: [dE/dx1s; dE/dx2s]
    hessE_local = zeros(6, 6); % Local Hessian
    potE = 0;                  % Potential Energy

    e1 = x1s - x2s; % Vector between points
    dist_sq = dot(e1, e1);
    I3 = eye(3);

    % Avoid division by zero or sqrt of zero if points coincide
    if dist_sq < 1e-24
        dist = 1e-12; % Use a very small number instead of zero
        e1_norm = [1;0;0]; % Assign arbitrary normalized direction
        dist_inv = 1 / dist;
        e1e1T_dist3 = zeros(3,3); % Avoid instability
    else
        dist = sqrt(dist_sq);
        dist_inv = 1 / dist;
        e1_norm = e1 * dist_inv; % Normalized direction vector
        e1e1T_dist3 = (e1 * e1') / (dist_sq * dist); % (e1*e1')/dist^3 term
    end
    ddist_dx1s = e1_norm; % d(dist)/d(x1s)
    d2dist_dx1s2 = dist_inv * I3 - e1e1T_dist3; % d^2(dist)/d(x1s)^2

    % Choose potential based on penetration state
    if state == 0 % NonPenetrated (Log-Exp Barrier)
        term1 = K1 * (h2 - dist);
        % Clamp exponent argument
        if term1 > 70; exp_term = exp(70); elseif term1 < -70; exp_term = exp(-70); else; exp_term = exp(term1); end

        one_plus_exp = 1 + exp_term;
        log_term = log(one_plus_exp);
        exp_over_one_plus_exp = exp_term / one_plus_exp;

        potE = (log_term / K1)^2;

        % dE/ddist = -(2/K1) * log_term * (K1 * exp_term / (1 + exp_term))
        dE_ddist = -2 * log_term * exp_over_one_plus_exp;

        % d2E/ddist2 = 2 * K1 * exp_over_one_plus_exp * (log_term * exp_term / one_plus_exp - exp_over_one_plus_exp)
        % Simplified from C++: d2E/ddist2 = 2*(exp/(1+exp))^2 + 2*log(1+exp)*(K1*exp/(1+exp)^2);
        d2E_ddist2 = 2 * exp_over_one_plus_exp^2 + 2 * log_term * (K1 * exp_term / (one_plus_exp^2));


    else % Penetrated (Quadratic Penalty)
        delta_dist = h2 - dist;
        potE = delta_dist^2;

        % dE/ddist = -2 * (h2 - dist)
        dE_ddist = -2 * delta_dist;
        % d2E/ddist2 = 2
        d2E_ddist2 = 2.0;
    end

    % Gradient: dE/dx = (dE/ddist) * (ddist/dx)
    gradE_local(1:3) = dE_ddist * ddist_dx1s;
    gradE_local(4:6) = -dE_ddist * ddist_dx1s; % ddist/dx2s = -ddist/dx1s

    % Hessian: d2E/dx2 = (d2E/ddist2)*(ddist/dx)*(ddist/dx)' + (dE/ddist)*d2dist/dx2
    hess_block = (d2E_ddist2 * (ddist_dx1s * ddist_dx1s')) + (dE_ddist * d2dist_dx1s2);
    hessE_local(1:3, 1:3) = hess_block;
    hessE_local(4:6, 4:6) = hess_block; % d2E/dx2s^2 = d2E/dx1s^2
    hessE_local(1:3, 4:6) = -hess_block; % d2E/dx1s dx2s = -d2E/dx1s^2
    hessE_local(4:6, 1:3) = -hess_block;


    % Map local to 12x1 gradient and 12x12 Hessian (P2P uses nodes 1 and 3 in 12x12 map)
    gradE = zeros(12, 1);
    hessE = zeros(12, 12);

    gradE(1:3) = gradE_local(1:3); % Node mapped to block 1 (idxA)
    gradE(7:9) = gradE_local(4:6); % Node mapped to block 3 (idxB)

    hessE(1:3, 1:3) = hessE_local(1:3, 1:3); % Block AA
    hessE(7:9, 7:9) = hessE_local(4:6, 4:6); % Block BB
    hessE(1:3, 7:9) = hessE_local(1:3, 4:6); % Block AB
    hessE(7:9, 1:3) = hessE_local(4:6, 1:3); % Block BA

    % Apply scaling
    gradE = gradE * contact_stiffness;
    hessE = hessE * contact_stiffness * scale; % Hessian scaling includes spatial derivative
    potE = potE * contact_stiffness * scale;   % Scale energy as well
end