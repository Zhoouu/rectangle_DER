function [dist, t, u, contactType, finalIndices] = computeSegmentSegmentDistance(p1, p2, p3, p4, idx1, idx2, idx3, idx4)
% Computes the minimum distance between two line segments (p1-p2 and p3-p4).
% Also determines the type of contact (Point-Point, Point-Edge, Edge-Edge)
% and the indices of the nodes involved in the closest feature pair.

    e1 = p2 - p1;
    e2 = p4 - p3;
    e12 = p3 - p1;

    D1 = dot(e1, e1);
    D2 = dot(e2, e2);
    R = dot(e1, e2);
    S1 = dot(e1, e12);
    S2 = dot(e2, e12);

    den = D1 * D2 - R^2;

    % Handle parallel or very short segments approximation
    if abs(den) < 1e-12 * D1 * D2 || D1 < 1e-12 || D2 < 1e-12
        % Simplified: Check endpoint-segment and endpoint-endpoint distances
        % This case needs more robust handling for production code,
        % checking all 4 endpoint-segment pairs and p1-p3, p1-p4, p2-p3, p2-p4 distances.
        % For now, just calculate based on t=0, u=0 assumption as a fallback.
        t = 0;
        u = -S2 / D2;
        u = max(0, min(1, u)); % Fix bound
        dist = norm(p1 - (p3 + u * e2));
        contactType = 1; % Defaulting to PointToEdge (idx1 vs edge idx2-(idx2+1))
        finalIndices = [idx1, idx2, idx1, idx2 + 1]; % Placeholder
        return; % Exit early for this simplified case
    end

    t = (S1 * D2 - S2 * R) / den;
    t_orig = t;
    t = max(0, min(1, t)); % Fix bound for t

    u = (t * R - S2) / D2;
    u_orig = u;
    u = max(0, min(1, u)); % Fix bound for u

    % Recompute t if u was clamped
    if abs(u - u_orig) > 1e-10
        t = (u * R + S1) / D1;
        t = max(0, min(1, t)); % Fix bound for t again
    end

    dist = norm((p1 + t * e1) - (p3 + u * e2));

    % Determine contact type and final indices
    is_t_endpoint = (abs(t) < 1e-9 || abs(t - 1) < 1e-9);
    is_u_endpoint = (abs(u) < 1e-9 || abs(u - 1) < 1e-9);

    if is_t_endpoint && is_u_endpoint % Point-Point
        contactType = 0; % PointToPoint Enum
        nodeA = idx1; if abs(t - 1) < 1e-9; nodeA = idx1 + 1; end
        nodeB = idx2; if abs(u - 1) < 1e-9; nodeB = idx2 + 1; end
        % For P2P, indices 3 and 4 are redundant/opposite node
        idx3_out = idx1 + 1; if abs(t - 1) < 1e-9; idx3_out = idx1; end
        idx4_out = idx2 + 1; if abs(u - 1) < 1e-9; idx4_out = idx2; end
        finalIndices = [nodeA, nodeB, idx3_out, idx4_out];

    elseif is_t_endpoint % Point (from edge 1) to Edge (edge 2)
        contactType = 1; % PointToEdge Enum
        nodeA = idx1; if abs(t - 1) < 1e-9; nodeA = idx1 + 1; end
        % For P2E, indices represent: Point A, Edge Start B, Point A, Edge End B
        finalIndices = [nodeA, idx2, nodeA, idx2 + 1];

    elseif is_u_endpoint % Point (from edge 2) to Edge (edge 1)
        contactType = 1; % PointToEdge Enum
        nodeB = idx2; if abs(u - 1) < 1e-9; nodeB = idx2 + 1; end
         % For P2E, indices represent: Edge Start A, Point B, Edge End A, Point B
        finalIndices = [idx1, nodeB, idx1 + 1, nodeB];

    else % Edge-Edge
        contactType = 2; % EdgeToEdge Enum
        finalIndices = [idx1, idx2, idx3, idx4];
    end

end