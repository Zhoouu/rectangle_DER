function simParams = defSimParams()
% This function defines the a struct contains the numerical parameters
%   Input:
%
%   Output:
%       simParams - the defined struct contains the numerical parameters 
%                   of the simulated system

simParams = struct();

simParams.totalTime = 9.0;

simParams.dt = 1e-2;

simParams.tol = 1e-6;

simParams.maximum_iter = 50;

simParams.Nsteps = round(simParams.totalTime/simParams.dt);

simParams.plotStep = 5;

% --- Added Contact Parameters ---
% Values below are examples, tune them based on C++ inputs & rod properties
simParams.col_limit = 0.0205;  % Broad phase distance limit (larger than 2*r0 + delta)
simParams.delta = 1e-6;     % Potential barrier thickness
simParams.k_scaler = 500;   % Contact stiffness scaler
simParams.mu = 0;         % Coefficient of friction (set to 0 to disable)
simParams.nu = 1e5;        % Friction velocity tolerance
% --- End Contact Parameters ---

end