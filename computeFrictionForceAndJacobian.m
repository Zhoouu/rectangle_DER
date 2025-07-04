function [fFriction, jFriction] = computeFrictionForceAndJacobian(rodParams, contact_gradient, indices, mu, K2, nu)
% Computes friction force and a simplified Jacobian.

    fFriction = zeros(12, 1);
    jFriction = zeros(12, 12); % Simplified Jacobian

    if mu == 0
        return; % No friction
    end

    dt = rodParams.dt;
    x = rodParams.x;
    x0 = rodParams.x0;

    idx1 = indices(1); idx2 = indices(2); idx3 = indices(3); idx4 = indices(4);

    % Get current and previous positions
    x1s = getVertex(x, idx1); x1e = getVertex(x, idx3);
    x2s = getVertex(x, idx2); x2e = getVertex(x, idx4);
    x1s0 = getVertex(x0, idx1); x1e0 = getVertex(x0, idx3);
    x2s0 = getVertex(x0, idx2); x2e0 = getVertex(x0, idx4);

    % Get contact forces on each node from the potential gradient
    % Note: contact_gradient is dE/dx, so force F = -dE/dx
    f1s = -contact_gradient(1:3); f1e = -contact_gradient(4:6);
    f2s = -contact_gradient(7:9); f2e = -contact_gradient(10:12);

    f1s_n = norm(f1s); f1e_n = norm(f1e);
    f2s_n = norm(f2s); f2e_n = norm(f2e);

    f1_tot = f1s + f1e;
    f2_tot = f2s + f2e; % Should be approx -f1_tot by Newton's 3rd law
    fn_mag = norm(f1_tot);

    if fn_mag < 1e-12 % Avoid division by zero if no contact force
        return;
    end

    contact_norm = f1_tot / fn_mag; % Normal direction based on force

    % Approximate interpolation factors (betas in C++)
    % These determine how velocity at nodes maps to velocity at contact point
    % Simple average for now, C++ uses force magnitudes which might be better
    beta11 = 0.5; beta12 = 0.5; % Assumes contact point is midpoint of edge 1 feature
    beta21 = 0.5; beta22 = 0.5; % Assumes contact point is midpoint of edge 2 feature
    % A better approximation might use the t, u parameters if available and reliable


    % Calculate velocities at nodes
    v1s = (x1s - x1s0) / dt; v1e = (x1e - x1e0) / dt;
    v2s = (x2s - x2s0) / dt; v2e = (x2e - x2e0) / dt;

    % Interpolate velocity at contact point on each body
    v1 = beta11 * v1s + beta12 * v1e;
    v2 = beta21 * v2s + beta22 * v2e;

    % Relative velocity
    v_rel = v1 - v2;

    % Tangential relative velocity
    tv_rel = v_rel - dot(v_rel, contact_norm) * contact_norm;
    tv_rel_n = norm(tv_rel);

    friction_type = 0; % 0:ZeroVel, 1:Sliding, 2:Sticking

    if tv_rel_n < 1e-12 % Nearly zero relative tangential velocity
        fFriction = zeros(12, 1); % No friction force
        friction_type = 0;
        % Jacobian is zero
    else
        tv_rel_u = tv_rel / tv_rel_n; % Unit tangential velocity direction

        if tv_rel_n > nu % Sliding regime
            gamma = 1.0;
            friction_type = 1;
        else % Sticking regime (smooth transition)
            gamma = (2.0 / (1.0 + exp(-K2 * tv_rel_n))) - 1.0;
            friction_type = 2;
        end

        % Friction force magnitude = mu * gamma * NormalForceMagnitude
        % We approximate NormalForceMagnitude by fn_mag
        ffr_val_vec = -mu * gamma * fn_mag * tv_rel_u; % Friction opposes relative tangential motion

        % Distribute friction force back to nodes (using same betas)
        fFriction(1:3)   = beta11 * ffr_val_vec;
        fFriction(4:6)   = beta12 * ffr_val_vec;
        fFriction(7:9)   = -beta21 * ffr_val_vec; % Acts on second body
        fFriction(10:12) = -beta22 * ffr_val_vec;

        % --- Simplified Friction Jacobian ---
        % The full analytic Jacobian is extremely complex.
        % Options:
        % 1. Zero Jacobian (simplest, may affect convergence)
        % 2. Identity matrix scaled by a factor (e.g., -mu * stiffness)
        % 3. Diagonal matrix based on damping analogy
        % We use Zero Jacobian here for simplicity. Replace if needed.
        jFriction = zeros(12, 12);
        % --- End Simplified Jacobian ---

    end
end