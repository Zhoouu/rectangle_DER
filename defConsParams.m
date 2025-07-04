function consParams = defConsParams(rodParams)
% This function defines the constrained parameters based on the physical
% struct's boundary conditions. 
%   Input:
%       rodParams - define a rod struct
%
%   Output:
%       consParams - the constrained parameters contains the constrained
%       index and unconstrained index of the discrete beam model

nv = rodParams.nv;
ne = rodParams.ne;

% Define fixed DOF
fixed_middle_rod_position=301:3*rodParams.nv;
fixed_middle_rod_theta=3*nv+99:rodParams.ndof;
fixed_start_end_position=1:12;
fixed_start_end_theta=3*nv+1:49:3*nv+50;
fixed_end_end_position=289:300;
fixed_end_end_theta=3*nv+49:49:3*nv+98;
% fixed_end_end_position=1177:1200;
% fixed_end_end_theta=3*nv+49:49:3*nv+392;
consInd = [fixed_start_end_position';fixed_end_end_position';fixed_start_end_theta';fixed_end_end_theta';fixed_middle_rod_position';fixed_middle_rod_theta'];

dummyInd = 1:rodParams.ndof;
dummyInd(consInd) = 0;
unconsInd = find(dummyInd>0);

consParams = struct();
consParams.unconsInd = unconsInd;
consParams.consInd = consInd;

end