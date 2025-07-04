function [Fc, Jc] = calculateContactForcesMATLAB(rodParams, simParams, contacts)
% Calculates total contact forces and Jacobians for all detected contacts.
% Uses updated potential helpers. Skips friction calculation if mu=0.

    ndof = rodParams.ndof;
    Fc = zeros(ndof, 1);
    Jc_rows = []; Jc_cols = []; Jc_vals = []; % For sparse Jacobian

    num_contacts = length(contacts);
    if num_contacts == 0
        Jc = sparse(ndof, ndof); return;
    end

    % Get parameters
    r0 = rodParams.r0;
    delta = simParams.delta;
    k_scaler = simParams.k_scaler;
    mu = simParams.mu; % Will be 0 based on defSimParams
    nu = simParams.nu;
    scale = 1 / r0;
    %scale = 1;
    h2 = 2 * r0 * scale;
    K1 = 250;
    % K1 = (15 * r0) / delta;
    K2 = 15 / nu; % (Not used if mu=0)

    % --- Contact Stiffness ---
    characteristic_elastic_force = rodParams.EA / rodParams.ne; % Approximation
    contact_stiffness = k_scaler * characteristic_elastic_force;

    for i = 1:num_contacts
        contact = contacts{i};
        indices = contact.indices;
        c_type = contact.type;
        c_state = contact.penetrationState;

        % Get relevant node positions (scaled for potential)
        x1s = getVertex(rodParams.x, indices(1)) * scale;
        x1e = getVertex(rodParams.x, indices(3)) * scale;
        x2s = getVertex(rodParams.x, indices(2)) * scale;
        x2e = getVertex(rodParams.x, indices(4)) * scale;

        % --- Calculate Potential Gradient and Hessian ---
        gradE_pot = zeros(12, 1);
        hessE_pot = zeros(12, 12); % Approximated Hessian

        if c_type == 0 % PointToPoint
             [gradE_pot, hessE_pot, ~] = computeContactPotentialP2P(x1s, x2s, K1, h2, c_state, contact_stiffness, scale);
        elseif c_type == 1 % PointToEdge
             % Determine which node is the point vs edge start/end
             idx_pt = -1; idx_e1 = -1; idx_e2 = -1;
             if indices(1) == indices(3); idx_pt=indices(1); idx_e1=indices(2); idx_e2=indices(4);
             elseif indices(2) == indices(4); idx_pt=indices(2); idx_e1=indices(1); idx_e2=indices(3);
             else; warning('Invalid P2E indices'); continue; end
             pt_coord = getVertex(rodParams.x, idx_pt) * scale;
             e1_coord = getVertex(rodParams.x, idx_e1) * scale;
             e2_coord = getVertex(rodParams.x, idx_e2) * scale;
             [gradE_mapped, hessE_mapped, ~] = computeContactPotentialP2E(e1_coord, e2_coord, pt_coord, K1, h2, c_state, contact_stiffness, scale);
             % Manually map back from [e1, e2, pt] order to [n1, n2, n3, n4] based on original indices
             map_from = {[1:3], [4:6], [7:9]}; % Blocks in P2E result [d/de1, d/de2, d/dpt]
             map_to_idx = {idx_e1, idx_e2, idx_pt}; % Corresponding node indices
             node_map_12x12 = {find(indices==map_to_idx{1}), find(indices==map_to_idx{2}), find(indices==map_to_idx{3})}; % Find 1,2,3,4 block for each node
             if length(node_map_12x12{1})>1; node_map_12x12{1}=node_map_12x12{1}(1); end % Handle P2P fallback case if needed
             if length(node_map_12x12{2})>1; node_map_12x12{2}=node_map_12x12{2}(1); end
             if length(node_map_12x12{3})>1; node_map_12x12{3}=node_map_12x12{3}(1); end

             for r=1:3
                 g_block_idx_r = (node_map_12x12{r}-1)*3 + (1:3);
                 gradE_pot(g_block_idx_r) = gradE_mapped(map_from{r});
                 for c=1:3
                     g_block_idx_c = (node_map_12x12{c}-1)*3 + (1:3);
                     hessE_pot(g_block_idx_r, g_block_idx_c) = hessE_mapped(map_from{r}, map_from{c});
                 end
             end

        elseif c_type == 2 % EdgeToEdge
             [gradE_pot, hessE_pot, ~] = computeContactPotentialE2E(x1s, x1e, x2s, x2e, K1, h2, c_state, contact_stiffness, scale);
        end

        % --- Friction Calculation (Skipped if mu=0) ---
        fFriction = zeros(12, 1);
        jFriction = zeros(12, 12);
        if mu > 0
            % Pass the potential hessian in case the friction jacobian needs it
            [fFriction, jFriction] = computeFrictionForceAndJacobian(rodParams, gradE_pot, hessE_pot, indices, mu, K2, nu);
        end
        % --- End Friction ---

        % --- Combine Force Gradient and Hessian ---
        total_grad = gradE_pot - fFriction; % Total Potential Gradient = Contact + Friction (if active)
        total_hess = hessE_pot - jFriction; % Total Hessian = Contact + Friction (if active)


        % --- Map local forces/Jacobian to Global System ---
        dof_map = @(node_idx) (3*(node_idx-1)+1):(3*(node_idx-1)+3);
        node_indices = indices;
        local_blocks = {[1:3], [4:6], [7:9], [10:12]};

        for row_block = 1:4
            node_r = node_indices(row_block);
            if node_r <= 0 || node_r > rodParams.nv; continue; end % Check valid node index
            global_row_dofs = dof_map(node_r);

            % Add force contribution (Force = -Gradient)
            Fc(global_row_dofs) = Fc(global_row_dofs) - total_grad(local_blocks{row_block});

            for col_block = 1:4
                node_c = node_indices(col_block);
                 if node_c <= 0 || node_c > rodParams.nv; continue; end % Check valid node index
                global_col_dofs = dof_map(node_c);

                % Add Jacobian contribution (Jacobian of Force = -Hessian of Potential)
                hess_block_3x3 = -total_hess(local_blocks{row_block}, local_blocks{col_block});

                [row_idx, col_idx] = ndgrid(global_row_dofs, global_col_dofs);
                Jc_rows = [Jc_rows; row_idx(:)];
                Jc_cols = [Jc_cols; col_idx(:)];
                Jc_vals = [Jc_vals; hess_block_3x3(:)];
            end
        end
    end % End loop over contacts

    Jc = sparse(Jc_rows, Jc_cols, Jc_vals, ndof, ndof);
end