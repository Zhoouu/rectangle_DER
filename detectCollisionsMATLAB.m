function contacts = detectCollisionsMATLAB(rodParams, simParams, edge)
% Finds potential and actual contacts between rod segments.

    % Parameters from structs
    r0 = rodParams.r0;
    x_current = rodParams.x;
    nv = rodParams.nv;
    ne = rodParams.ne;
    delta = simParams.delta;
    col_limit = simParams.col_limit;

    % Calculate derived parameters
    scale = 1 / r0;
    % scale = 1;
    contact_limit_sq = (scale * (2 * r0 + delta))^2;     % Squared for efficiency
    candidate_limit_sq = (scale * (2 * r0 + col_limit))^2; % Squared for efficiency
    numerical_limit = scale * (2 * r0 - delta);
    surface_limit_sq = (scale * 2 * r0)^2;             % Squared for efficiency

    ignore_adjacent = 3; % As per C++ code

    candidate_list = {}; % Cell array to store potential contacts
    num_candidates = 0;

    % --- Phase 1: Broad Phase - Construct Candidate Set ---
    for i = 1:ne
        for j = (i + ignore_adjacent + 1):ne
            % Get edge endpoints (scaled)
            idx1 = edge(i, 1); idx1_p1 = edge(i, 2);
            idx2 = edge(j, 1); idx2_p1 = edge(j, 2);

            p1 = getVertex(x_current, idx1) * scale;
            p2 = getVertex(x_current, idx1_p1) * scale;
            p3 = getVertex(x_current, idx2) * scale;
            p4 = getVertex(x_current, idx2_p1) * scale;

            % Simple bounding box check first (optional but can speed up)
            min1 = min(p1, p2); max1 = max(p1, p2);
            min2 = min(p3, p4); max2 = max(p3, p4);
            box_dist_sq = sum(max(0, max(min1 - max2, min2 - max1)).^2);
            if box_dist_sq > candidate_limit_sq
                continue; % Boxes too far apart
            end

            % Compute rough distance (more precise than bbox)
            [dist, ~, ~, ~, ~] = computeSegmentSegmentDistance(p1, p2, p3, p4, idx1, idx2, idx1_p1, idx2_p1);
            dist_sq = dist^2;

            % --- Early Exit Check (Penetration during candidate search) ---
            % C++ code had an option `ignore_escape`. We'll assume it's false.
            % if dist_sq < surface_limit_sq
            %     warning('Rod penetration detected during candidate search between edges %d and %d!', i, j);
            %     % Handle this case - maybe stop simulation or log error
            %     % For now, we'll let it proceed to narrow phase
            % end
            % --- End Early Exit Check ---


            if dist_sq < candidate_limit_sq
                num_candidates = num_candidates + 1;
                candidate_list{num_candidates} = [i, j]; % Store original edge indices
            end
        end
    end

    % --- Phase 2: Narrow Phase - Refine Candidates ---
    contacts = {}; % Cell array for final contacts
    num_collisions = 0;

    for k = 1:num_candidates
        edge_idx1 = candidate_list{k}(1);
        edge_idx2 = candidate_list{k}(2);

        idx1 = edge(edge_idx1, 1); idx1_p1 = edge(edge_idx1, 2);
        idx2 = edge(edge_idx2, 1); idx2_p1 = edge(edge_idx2, 2);

        p1 = getVertex(x_current, idx1) * scale;
        p2 = getVertex(x_current, idx1_p1) * scale;
        p3 = getVertex(x_current, idx2) * scale;
        p4 = getVertex(x_current, idx2_p1) * scale;

        [dist, t, u, contactType, finalIndices] = computeSegmentSegmentDistance(p1, p2, p3, p4, idx1, idx2, idx1_p1, idx2_p1);
        
        % --- MODIFICATION START: Only process E2E contacts ---
        if contactType ~= 2 % Skip if not EdgeToEdge (type 2)
            continue; % Go to the next candidate
        end
        % --- MODIFICATION END ---
        
        dist_sq = dist^2;

        if dist_sq < contact_limit_sq
            num_collisions = num_collisions + 1;
            contact_info = struct();
            contact_info.indices = finalIndices; % [idxA, idxB, idxC, idxD] from helper
            contact_info.type = contactType;     % 0:P2P, 1:P2E, 2:E2E
            contact_info.distance = dist / scale; % Unscaled distance
            contact_info.params = [t, u];

            if dist > numerical_limit
                contact_info.penetrationState = 0; % NonPenetrated
            else
                contact_info.penetrationState = 1; % Penetrated
            end
            contacts{num_collisions} = contact_info;
        end
    end
    % fprintf('Detected %d contacts.\n', num_collisions);
end