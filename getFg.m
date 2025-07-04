function [Fg, Jg] = getFg(rodParams)
% This function computes the gravitational force and jacobian of the simulated system.
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%
% OUTPUTS: Fg - gravitational forces (ndof x 1)
%          Jg - gravitational jacobian (ndof x ndof)
global global_x0;
garr = zeros(rodParams.ndof, 1);
for c = 1:rodParams.nv
    garr( 3 * (c-1) + 1 : 3 * (c-1) + 3) = [0;global_x0(3 * (c-1) + 2);global_x0(3 * (c-1) + 3)];
end
garr=garr/norm(garr)*norm(rodParams.g);
Fg = rodParams.m .* garr;
Jg = zeros(rodParams.ndof, rodParams.ndof);

end
