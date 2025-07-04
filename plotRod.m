function plotRod(x, nv, edge, edgeColors)
% This function visualizes the rod by plotting its edges as lines.
%
%   Input:
%       x    - Global DOF vector (ndof x 1)
%       nv   - number of vertices
%       edge - Edge connectivity list

useDefaultColor = false;
if nargin < 4
    useDefaultColor = true;
end

h1 = figure(1);

clf()
hold on

% Loop through each edge defined in the connectivity list
for i = 1:size(edge, 1)
    % Get the indices of the two nodes forming the edge
    node_idx1 = edge(i, 1);
    node_idx2 = edge(i, 2);

    % Extract the 3D coordinates for each node using the helper function
    pos1 = getVertex(x, node_idx1);
    pos2 = getVertex(x, node_idx2);

    % 确定当前边的颜色
    if useDefaultColor
        lineColor = 'k'; % 使用默认黑色
    else
        lineColor = edgeColors{i}; % 使用指定的颜色
    end

    % Plot a line segment (the edge) between the two nodes
    plot3([pos1(1), pos2(1)], [pos1(2), pos2(2)], [pos1(3), pos2(3)], 'Color', lineColor, 'LineWidth', 2);
end

% --- Keep existing axis and view settings ---
axis([-1 1 -1 1 -2 1]/2);
view([1 1 1]);
camup([1 0 0]);
xlim([-0.15 0.65]);
zlim([-0.5 0.5]);

hold off
drawnow

end