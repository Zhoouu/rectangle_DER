function plotRod(x, nv)
% This function visualizes the rod by plotting its nodal positions.
%
%   Input:
%       x - Global DOF vector (ndof x 1)

x1 = x(1:3:3*nv-2);
x2 = x(2:3:3*nv-1);
x3 = x(3:3:3*nv);

h1 = figure(1);

clf()
plot3(x1,x2, x3, '-');
hold on
% axis([-1 1 -1 1 -1.5 0.5]);
% view(90,45);
axis([-1 1 -1 1 -2 1]/2);
view([1 1 1]);
camup([1 0 0]);
xlim([-0.15 0.65]);        % 只改 x 轴
%ylim([-4  8]);         % 只改 y 轴
zlim([-0.5 0.5]);         % 只改 z 轴

drawnow
