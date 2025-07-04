function [xUncons, iter, normf] = objfun(rodParams, simParams, consParams, sElement, bElement,edge)
% This function solves for the equilibrium/dynamic state of a simulated system
%        using Newton's method on unconstrained degrees of freedom.
%
%   Input:
%       rodParams  - Struct containing rod material and state parameters
%       simParams  - Struct with numerical parameters
%       consParams - Struct with constraints (unconstrained DOF indices)
%       sElement   - Stretching element definitions
%       bElement   - Bending element definitions
%       edge       - Edge connectivity list (needed for collision detection)
%
%   Output:
%       xUncons    - Updated positions of unconstrained DOFs after convergence

% Numerical parameter
maximum_iter = simParams.maximum_iter;
tol = simParams.tol;
dt = simParams.dt;

% Free DOF
unconsInd = consParams.unconsInd;
ndof_uncons = length(unconsInd); % Number of unconstrained DOFs

% Mass matrix
mUncons = rodParams.m(unconsInd);
mMat = diag(mUncons);

% Damping
viscosity = rodParams.viscosity;

% DOF vector
xUncons = rodParams.x(unconsInd);
uUncons = rodParams.u(unconsInd);

% previous position
x0Uncons = rodParams.x0(consParams.unconsInd);

% number of iterations
iter = 0;

% norm of function value
normf = tol*10;
normf0 = -1; % Store initial norm

alpha = 1; % Damping factor for Newton step (optional)

% Newton's method
while (normf > tol)
    if iter > 0 && normf <= tol * normf0 % Relative tolerance check
        fprintf('Converged with relative tolerance.\n');
        break;
    end

    rodParams.x(unconsInd) = xUncons;

    % Update stretching element
    for i = 1:rodParams.ne
        sElement(i).nodePos_1 = getVertex(rodParams.x, sElement(i).nodeIndex(1));
        sElement(i).nodePos_2 = getVertex(rodParams.x, sElement(i).nodeIndex(2));
        % Update tangent if needed by other calculations
        edge_vec = sElement(i).nodePos_2 - sElement(i).nodePos_1;
        sElement(i).t = edge_vec / norm(edge_vec);
    end

    % Update bending element
    for i =1:rodParams.nb
        bElement(i).nodePos_1 = getVertex(rodParams.x, bElement(i).nodeIndex(1));
        bElement(i).nodePos_2 = getVertex(rodParams.x, bElement(i).nodeIndex(2));
        bElement(i).nodePos_3 = getVertex(rodParams.x, bElement(i).nodeIndex(3));

        if (bElement(i).directSign_1 > 0)
            bElement(i).theta_1 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        else
            bElement(i).theta_1 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        end

        if (bElement(i).directSign_2 > 0)
            bElement(i).theta_2 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        else
            bElement(i).theta_2 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        end

        % update frame (tangents first)
        t1_new = (bElement(i).nodePos_2 - bElement(i).nodePos_1);
        t2_new = (bElement(i).nodePos_3 - bElement(i).nodePos_2);
        norm_t1 = norm(t1_new); norm_t2 = norm(t2_new);
        if norm_t1 < 1e-12 || norm_t2 < 1e-12
            warning('Near zero edge length in bending element %d', i);
            % Keep old frames if edge length is near zero
            bElement(i).t_1 = bElement(i).t_1_old;
            bElement(i).t_2 = bElement(i).t_2_old;
            bElement(i).d_11 = bElement(i).d_11_old;
            bElement(i).d_12 = bElement(i).d_12_old;
            bElement(i).d_21 = bElement(i).d_21_old;
            bElement(i).d_22 = bElement(i).d_22_old;
            bElement(i).refTwist = bElement(i).refTwist_old;
        else
            bElement(i).t_1 = t1_new / norm_t1;
            bElement(i).t_2 = t2_new / norm_t2;

            % Parallel transport directors
            bElement(i).d_11 = parallel_transport(bElement(i).d_11_old, bElement(i).t_1_old, bElement(i).t_1);
            bElement(i).d_11 = bElement(i).d_11 / norm(bElement(i).d_11); % Re-normalize
            bElement(i).d_12 = cross(bElement(i).t_1, bElement(i).d_11);
            bElement(i).d_12 = bElement(i).d_12 / norm(bElement(i).d_12); % Not strictly necessary if d_11 and t_1 are ortho-normal

            bElement(i).d_21 = parallel_transport(bElement(i).d_21_old, bElement(i).t_2_old, bElement(i).t_2);
            bElement(i).d_21 = bElement(i).d_21 / norm(bElement(i).d_21); % Re-normalize
            bElement(i).d_22 = cross(bElement(i).t_2, bElement(i).d_21);
            bElement(i).d_22 = bElement(i).d_22 / norm(bElement(i).d_22);

            % Compute reference twist update (important!)
            u_temp = parallel_transport(bElement(i).d_11, bElement(i).t_1, bElement(i).t_2); % Transport d11 along segment
            u_temp = rotateAxisAngle(u_temp, bElement(i).t_2, bElement(i).refTwist_old); % Apply old twist
            deltaAngle = signedAngle(u_temp, bElement(i).d_21, bElement(i).t_2); % Angle between transported old d1 and new d1
            bElement(i).refTwist = bElement(i).refTwist_old + deltaAngle; % Update twist angle

        end

        % build material frame
        cs = cos( bElement(i).theta_1 );
        ss = sin( bElement(i).theta_1 );
        bElement(i).m_11 =   cs * bElement(i).d_11 + ss * bElement(i).d_12;
        bElement(i).m_12 = - ss * bElement(i).d_11 + cs * bElement(i).d_12;

        cs = cos( bElement(i).theta_2 );
        ss = sin( bElement(i).theta_2 );
        bElement(i).m_21 =   cs * bElement(i).d_21 + ss * bElement(i).d_22;
        bElement(i).m_22 = - ss * bElement(i).d_21 + cs * bElement(i).d_22;
    end

    % --- Start Contact Calculation ---
    contacts = detectCollisionsMATLAB(rodParams, simParams, edge);
    [Fc, Jc] = calculateContactForcesMATLAB(rodParams, simParams, contacts);
    % --- End Contact Calculation ---

    % Get forces
    [Fg, Jg] = getFg(rodParams);
    [Fs, Js] = getFs(rodParams, sElement);
    [Fb, Jb] = getFb(rodParams, bElement);
    [Ft, Jt] = getFt(rodParams, bElement);

    % --- Assemble Total Forces and Jacobian ---
    Forces_internal = Fg + Fs + Fb + Ft; % Sum of internal/gravity forces
    Forces_contact = Fc;                % Contact forces
    TotalForces = Forces_internal + Forces_contact; % Total external/internal/contact forces

    JForces_internal = Jg + Js + Jb + Jt; % Jacobian of internal/gravity forces
    JContact = Jc;                        % Jacobian of contact forces
    TotalJacobian = JForces_internal + JContact; % Total Jacobian

    % Extract unconstrained parts
    ForcesUncons = TotalForces(unconsInd);
    JacobianUncons = TotalJacobian(unconsInd, unconsInd);
    % --- End Assembly ---

    % --- Equation of Motion Residual (f = M*a + C*v - F_ext = 0) ---
    if (rodParams.ifStatic == 1)
        f = - ForcesUncons; % Static equilibrium: Sum of forces = 0
        J = - JacobianUncons; % Jacobian of the residual w.r.t xUncons
    else
        % Dynamic implicit step (Backward Euler / Implicit Midpoint style)
        % Residual f = M*(x - x0 - u0*dt)/dt^2 + C*(x-x0)/dt - F_total(x) = 0
        inertial_term_f = mUncons .* ( (xUncons - x0Uncons)/dt^2 - uUncons/dt );
        damping_term_f = viscosity * mUncons .* (xUncons - x0Uncons)/dt;
        f = inertial_term_f + damping_term_f - ForcesUncons;

        % Jacobian J = d(f)/d(xUncons)
        inertial_term_J = mMat/dt^2;
        damping_term_J = viscosity * mMat/dt;
        J = inertial_term_J + damping_term_J - JacobianUncons;
    end
    % --- End Equation of Motion ---

    % Calculate norm of the residual
    normfNew = norm(f);
    if iter == 0
        normf0 = normfNew; % Store initial norm
        if normf0 < tol
            fprintf('Initial guess already converged.\n');
            break; % Already converged
        end
    end

    fprintf('Iter=%d, error=%e, alpha=%f\n', iter, normfNew, alpha);

    % Newton's update step (solve J*dx = -f)
    if rcond(full(J)) < 1e-14 % Check condition number before solving
        warning('Jacobian is ill-conditioned (rcond = %e). Newton step might be inaccurate.', rcond(full(J)));
        % Consider alternatives: different solver, regularization, smaller dt, or stop.
        % For now, we attempt the solve anyway.
        dx = -J\f;
    else
        dx = -J\f;
    end


    % --- Line Search / Damping (Optional but often helpful) ---
    % Simple backtracking line search can be added here if needed
    % Or just simple damping:
    % alpha = min(1.0, max(0.1, alpha * 0.95)); % Simple adaptive damping
    %alpha = 1.0; % No damping for now
    % --- End Line Search ---

    % Newton's update
    xUncons = xUncons + alpha * dx; % Update the guess

    % Update iteration count and old norm
    iter = iter + 1;
    normf = normfNew;

    if (iter > maximum_iter)
        warning('Newton method did not converge within %d iterations. Norm = %e', maximum_iter, normf);
        break; % Exit loop if max iterations reached
    end
end

rodParams.x(unconsInd) = xUncons;

end

