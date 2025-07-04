function [gradE, hessE, potE] = computeContactPotentialE2E(x1s, x1e, x2s, x2e, K1, h2, state, contact_stiffness, scale)
% Calculates E2E potential gradient, Hessian (Approximated Analytic), and Energy.
% Implements Analytic Gradient. Hessian uses approximated analytic form.
% Includes explicit reshaping and checks for vector dimensions.

    gradE = zeros(12, 1); % Global gradient [dE/dx1s; dE/dx1e; dE/dx2s; dE/dx2e]
    hessE = zeros(12, 12);% Global Hessian (Approximated Analytic)
    potE = 0;             % Potential Energy

    % --- Calculate Analytic Gradient, Energy and required derivatives ---
    [gradE, potE, dE_ddist_abs, dist, ddist_dx_vec, d2E_ddist2] = computeE2E_Gradient_Energy_Derivs(x1s, x1e, x2s, x2e, K1, h2, state);


    % --- Hessian Approximation ---
    % Full Hessian: d2E/dx2 = (d2E/d|dist|^2)*(d|dist|/dx)*(d|dist|/dx)' + (dE/d|dist|)*d2|dist|/dx2
    % Approximation: Neglect the second term involving d2|dist|/dx2 (complex).
    if dist > 1e-12 % Avoid issues if distance is zero
        % Ensure ddist_dx_vec is column vector before outer product
        if size(ddist_dx_vec, 2) ~= 1
            warning('E2E Hessian: ddist_dx_vec was not a column vector. Reshaping.');
            ddist_dx_vec = ddist_dx_vec(:); % Force into column vector
        end
         if size(ddist_dx_vec, 1) ~= 12
             error('E2E Hessian: ddist_dx_vec does not have 12 elements.');
         end
        hessE = d2E_ddist2 * (ddist_dx_vec * ddist_dx_vec');
    else
        hessE = zeros(12,12); % Hessian is zero if distance is zero
    end
    % Symmetrize (should already be symmetric by construction here)
    % hessE = 0.5 * (hessE + hessE');

    % --- Scale ---
    gradE = gradE * contact_stiffness;
    hessE = hessE * contact_stiffness * scale; % Hessian scaling
    potE = potE * contact_stiffness * scale;   % Scale energy
end


% --- Helper function for E2E Gradient/Energy/Derivatives ---
function [gradE, potE, dE_ddist_abs, dist, ddist_dx_vec, d2E_ddist2] = computeE2E_Gradient_Energy_Derivs(x1s, x1e, x2s, x2e, K1, h2, state)
    % Calculates Energy, Gradient, dE/d|dist|, distance, d|dist|/dx, and d2E/d|dist|^2
    gradE = zeros(12, 1);
    potE = 0;
    dE_ddist_abs = 0; % Derivative w.r.t absolute distance |dist|
    dist = 0;
    ddist_dx_vec = zeros(12,1);
    d2E_ddist2 = 0; % Second derivative of potential w.r.t absolute distance
    I3 = eye(3);

    e1 = x1e - x1s;
    e2 = x2e - x2s;
    e13 = x2s - x1s;

    n = cross(e1, e2);
    n_norm_sq = dot(n, n);

    % Initialize derivatives to zero vectors
    ddist_signed_dx1s = zeros(3,1); ddist_signed_dx1e = zeros(3,1);
    ddist_signed_dx2s = zeros(3,1); ddist_signed_dx2e = zeros(3,1);
    sgn_dot = 1; % Default sign

    if n_norm_sq < 1e-24 % Edges are parallel
       % warning('E2E Inner: Parallel edges. Using P2P midpoint approx.');
       mid1 = (x1s + x1e) / 2;
       mid2 = (x2s + x2e) / 2;
       dist_vec = mid1 - mid2;
       dist_sq_p2p = dot(dist_vec, dist_vec);
       if dist_sq_p2p < 1e-24; dist=1e-12; else; dist=sqrt(dist_sq_p2p); end

       if state == 0 % NonPenetrated
            term1 = K1 * (h2 - dist);
            if term1 > 70; exp_term = exp(70); elseif term1 < -70; exp_term = exp(-70); else; exp_term = exp(term1); end
            one_plus_exp = 1 + exp_term; log_term = log(one_plus_exp); exp_over_one_plus_exp = exp_term / one_plus_exp;
            potE = (log_term / K1)^2;
            dE_ddist_abs = -2 * log_term * exp_over_one_plus_exp;
            d2E_ddist2 = 2 * exp_over_one_plus_exp^2 + 2 * log_term * (K1 * exp_term / (one_plus_exp^2));
       else % Penetrated
            potE = (h2 - dist)^2;
            dE_ddist_abs = -2 * (h2 - dist);
            d2E_ddist2 = 2.0;
       end
       e_norm = dist_vec/max(dist,1e-12);
       ddist_dmid1 = e_norm; ddist_dmid2 = -e_norm;
       ddist_signed_dx1s = 0.5 * ddist_dmid1; ddist_signed_dx1e = 0.5 * ddist_dmid1;
       ddist_signed_dx2s = 0.5 * ddist_dmid2; ddist_signed_dx2e = 0.5 * ddist_dmid2;
       sgn_dot = 1; % Treat as non-signed case for gradient calculation below
    else
        n_norm = sqrt(n_norm_sq);
        n_unit = n / n_norm;
        n_norm_inv = 1 / n_norm;
        dot_n_e13 = dot(n, e13);
        dist = abs(dot_n_e13) * n_norm_inv; % Absolute distance

        sgn_dot = sign(dot_n_e13);
        if sgn_dot == 0; sgn_dot = 1; end % Handle exact zero case

        % Compute Potential, dE/d|dist|, d2E/d|dist|^2 based on state
        if state == 0 % NonPenetrated
            term1 = K1 * (h2 - dist);
            if term1 > 70; exp_term = exp(70); elseif term1 < -70; exp_term = exp(-70); else; exp_term = exp(term1); end
            one_plus_exp = 1 + exp_term; log_term = log(one_plus_exp); exp_over_one_plus_exp = exp_term / one_plus_exp;
            potE = (log_term / K1)^2;
            dE_ddist_abs = -2 * log_term * exp_over_one_plus_exp;
            d2E_ddist2 = 2 * exp_over_one_plus_exp^2 + 2 * log_term * (K1 * exp_term / (one_plus_exp^2));
        else % Penetrated
            delta_dist = h2 - dist;
            potE = delta_dist^2;
            dE_ddist_abs = -2 * delta_dist;
            d2E_ddist2 = 2.0;
        end

        % Analytic Gradient of signed distance P/Q = dot(n,e13)/norm(n)
        e1_x = crossMat(e1); e2_x = crossMat(e2);
        dn_dx1s = crossMat(e2); dn_dx1e = -crossMat(e2); dn_dx2s = -crossMat(e1); dn_dx2e = crossMat(e1);

        dQ_dx1s = (n_unit' * dn_dx1s)'; dQ_dx1e = (n_unit' * dn_dx1e)';
        dQ_dx2s = (n_unit' * dn_dx2s)'; dQ_dx2e = (n_unit' * dn_dx2e)';

        de13_dx1s_mat = -I3; de13_dx1e_mat = zeros(3,3); de13_dx2s_mat = I3; de13_dx2e_mat = zeros(3,3);

        dP_dx1s = (dn_dx1s' * e13) + (n' * de13_dx1s_mat)';
        dP_dx1e = (dn_dx1e' * e13) + (n' * de13_dx1e_mat)';
        dP_dx2s = (dn_dx2s' * e13) + (n' * de13_dx2s_mat)';
        dP_dx2e = (dn_dx2e' * e13) + (n' * de13_dx2e_mat)';

        % Force results to be 3x1 column vectors
        dP_dx1s = dP_dx1s(:); dP_dx1e = dP_dx1e(:); dP_dx2s = dP_dx2s(:); dP_dx2e = dP_dx2e(:);
        dQ_dx1s = dQ_dx1s(:); dQ_dx1e = dQ_dx1e(:); dQ_dx2s = dQ_dx2s(:); dQ_dx2e = dQ_dx2e(:);

        % Gradient of signed distance P/Q
        ddist_signed_dx1s = (dP_dx1s * n_norm - dot_n_e13 * dQ_dx1s) / n_norm_sq;
        ddist_signed_dx1e = (dP_dx1e * n_norm - dot_n_e13 * dQ_dx1e) / n_norm_sq;
        ddist_signed_dx2s = (dP_dx2s * n_norm - dot_n_e13 * dQ_dx2s) / n_norm_sq;
        ddist_signed_dx2e = (dP_dx2e * n_norm - dot_n_e13 * dQ_dx2e) / n_norm_sq;
    end % End if not parallel

    % Final Gradient of Potential: dE/dx = (dE/d|dist|) * d|dist|/dx
    % d|dist|/dx = d(abs(P/Q))/dx = sgn(P/Q) * d(P/Q)/dx = sgn_dot * d(signed_dist)/dx
    ddist_dx1s = sgn_dot * ddist_signed_dx1s;
    ddist_dx1e = sgn_dot * ddist_signed_dx1e;
    ddist_dx2s = sgn_dot * ddist_signed_dx2s;
    ddist_dx2e = sgn_dot * ddist_signed_dx2e;

    % Assemble 12x1 gradient of absolute distance
    ddist_dx_vec = [ddist_dx1s; ddist_dx1e; ddist_dx2s; ddist_dx2e];

    % --- Final Gradient of Potential Energy ---
    gradE = dE_ddist_abs * ddist_dx_vec;

    % --- Final Sanity Check for Dimensions before returning ---
    if ~isequal(size(gradE),[12,1])
        error('E2E Final Gradient calculation resulted in wrong dimensions: %s', mat2str(size(gradE)));
    end
     if ~isequal(size(ddist_dx_vec),[12,1])
        error('E2E Final ddist_dx_vec calculation resulted in wrong dimensions: %s', mat2str(size(ddist_dx_vec)));
    end

end