function [gradE, hessE, potE] = computeContactPotentialP2E(x1s, x1e, x2s, K1, h2, state, contact_stiffness, scale)
% Calculates P2E potential gradient, Hessian (Approximated Analytic), and Energy.
% Hessian Approximation: Hess ~ (d2E/ddist2) * (ddist/dx)(ddist/dx)'

    gradE_local = zeros(9, 1); % Local gradient: [dE/dx1s; dE/dx1e; dE/dx2s]
    hessE_local = zeros(9, 9); % Local Hessian (Approximated Analytic)
    potE = 0;                  % Potential Energy

    % --- Calculate Analytic Gradient, Energy, and required derivatives ---
    [gradE_local, potE, dE_ddist, dist, ddist_dx_vec, d2E_ddist2] = computeP2E_Gradient_Energy_Derivs(x1s, x1e, x2s, K1, h2, state);

    % --- Hessian Approximation ---
    % Full Hessian: d2E/dx2 = (d2E/ddist2)*(ddist/dx)*(ddist/dx)' + (dE/ddist)*d2dist/dx2
    % Approximation: Neglect the second term involving d2dist/dx2 (complex).
    if dist > 1e-12 % Avoid issues if distance is zero
        hessE_local = d2E_ddist2 * (ddist_dx_vec * ddist_dx_vec');
    else
        hessE_local = zeros(9,9); % Hessian is zero if distance is zero
    end
    % Symmetrize (should already be symmetric by construction here)
    % hessE_local = 0.5 * (hessE_local + hessE_local');


    % Map local 9x1 grad / 9x9 hess to 12x1 / 12x12 output
    gradE = zeros(12, 1);
    hessE = zeros(12, 12);
    % Mapping assumes x1s->Block 1, x1e->Block 2, x2s->Block 3
    indices_map_9 = {[1:3], [4:6], [7:9]}; % Local 9x9 blocks
    indices_map_12 = {[1:3], [4:6], [7:9]}; % Corresponding global 12x12 blocks

    for r = 1:3
        gradE(indices_map_12{r}) = gradE_local(indices_map_9{r});
        for c = 1:3
            hessE(indices_map_12{r}, indices_map_12{c}) = hessE_local(indices_map_9{r}, indices_map_9{c});
        end
    end

    % Apply scaling
    gradE = gradE * contact_stiffness;
    hessE = hessE * contact_stiffness * scale; % Hessian scaling
    potE = potE * contact_stiffness * scale;   % Scale energy
end

% --- Helper function for P2E Gradient/Energy/Derivatives ---
function [gradE_local, potE, dE_ddist, dist, ddist_dx_vec, d2E_ddist2] = computeP2E_Gradient_Energy_Derivs(x1s, x1e, x2s, K1, h2, state)
    % Calculates Energy, Gradient, dE/ddist, distance, d(dist)/dx, and d2E/ddist2
    gradE_local = zeros(9, 1);
    potE = 0;
    dE_ddist = 0;
    dist = 0;
    ddist_dx_vec = zeros(9,1);
    d2E_ddist2 = 0; % Second derivative of potential w.r.t distance

    e1  = x1e - x1s;  % Edge vector
    e13 = x1s - x2s;
    e23 = x1e - x2s;

    edgeLen_sq = dot(e1, e1);
     if edgeLen_sq < 1e-24
        % Degenerate edge - treat as P2P between x1s and x2s
        % warning('P2E Inner: Degenerate edge. Using P2P.');
        dist_vec = x1s - x2s;
        dist_sq_p2p = dot(dist_vec, dist_vec);
        if dist_sq_p2p < 1e-24; dist = 1e-12; else; dist = sqrt(dist_sq_p2p); end

        if state == 0 % NonPenetrated
            term1 = K1 * (h2 - dist);
            if term1 > 70; exp_term = exp(70); elseif term1 < -70; exp_term = exp(-70); else; exp_term = exp(term1); end
            one_plus_exp = 1 + exp_term; log_term = log(one_plus_exp); exp_over_one_plus_exp = exp_term / one_plus_exp;
            potE = (log_term / K1)^2;
            dE_ddist = -2 * log_term * exp_over_one_plus_exp;
            d2E_ddist2 = 2 * exp_over_one_plus_exp^2 + 2 * log_term * (K1 * exp_term / (one_plus_exp^2));
        else % Penetrated
            potE = (h2 - dist)^2;
            dE_ddist = -2 * (h2 - dist);
            d2E_ddist2 = 2.0;
        end
        e1_norm = dist_vec/max(dist,1e-12);
        ddist_dx1s_p2p = e1_norm;
        ddist_dx2s_p2p = -e1_norm;
        gradE_local(1:3) = dE_ddist * ddist_dx1s_p2p; % dE/dx1s
        gradE_local(7:9) = dE_ddist * ddist_dx2s_p2p; % dE/dx2s
        ddist_dx_vec(1:3) = ddist_dx1s_p2p;
        ddist_dx_vec(7:9) = ddist_dx2s_p2p;
        return
    end
    edgeLen = sqrt(edgeLen_sq);
    edgeLen_inv = 1 / edgeLen;

    v = cross(e13, e23);
    v_norm_sq = dot(v, v);
    if v_norm_sq < 1e-24
        dist = 0;
        v_norm = 1e-12; v_unit = [1;0;0];
    else
        v_norm = sqrt(v_norm_sq); v_unit = v / v_norm;
        dist = v_norm * edgeLen_inv;
    end

    % Compute Potential, dE/ddist, d2E/ddist2 based on state
    if state == 0 % NonPenetrated
        term1 = K1 * (h2 - dist);
        if term1 > 70; exp_term = exp(70); elseif term1 < -70; exp_term = exp(-70); else; exp_term = exp(term1); end
        one_plus_exp = 1 + exp_term; log_term = log(one_plus_exp); exp_over_one_plus_exp = exp_term / one_plus_exp;
        potE = (log_term / K1)^2;
        dE_ddist = -2 * log_term * exp_over_one_plus_exp;
        d2E_ddist2 = 2 * exp_over_one_plus_exp^2 + 2 * log_term * (K1 * exp_term / (one_plus_exp^2));
    else % Penetrated
        delta_dist = h2 - dist;
        potE = delta_dist^2;
        dE_ddist = -2 * delta_dist;
        d2E_ddist2 = 2.0;
    end

    % Analytic Gradient of Distance d(dist)/dx
    e1_norm = e1 * edgeLen_inv;
    dv_dx1s = -crossMat(e23); dv_dx1e = crossMat(e13); dv_dx2s = crossMat(e23) - crossMat(e13);
    dN_dx1s = (v_unit' * dv_dx1s)'; dN_dx1e = (v_unit' * dv_dx1e)'; dN_dx2s = (v_unit' * dv_dx2s)';
    dD_dx1s = -e1_norm; dD_dx1e = e1_norm; dD_dx2s = zeros(3,1);

    if edgeLen_sq < 1e-24 % Should be caught above, but re-check
        ddist_dx1s = zeros(3,1); ddist_dx1e = zeros(3,1); ddist_dx2s = zeros(3,1);
    else
        ddist_dx1s = (dN_dx1s * edgeLen - v_norm * dD_dx1s) / edgeLen_sq;
        ddist_dx1e = (dN_dx1e * edgeLen - v_norm * dD_dx1e) / edgeLen_sq;
        ddist_dx2s = (dN_dx2s * edgeLen - v_norm * dD_dx2s) / edgeLen_sq;
    end
    ddist_dx_vec = [ddist_dx1s; ddist_dx1e; ddist_dx2s]; % Store 9x1 gradient of distance

    % Gradient dE/dx = (dE/ddist) * (ddist/dx)
    gradE_local = dE_ddist * ddist_dx_vec;
end