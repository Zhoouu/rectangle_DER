function [Fs, Js] = getFs(rodParams, sElement)
% This function computes the stretching force and jacobian of the simulated system. 
% INPUTS: rodParams - the defined rod struct contains the physical and
%                     numerical parameters of the simulated system
%         sElement - the stretching element list (ne x 2)
%
% OUTPUTS: Fs - stretching forces (ndof x 1)
%          Js - stretching jacobian (ndof x ndof)

Fs = zeros(rodParams.ndof, 1);
Js = zeros(rodParams.ndof, rodParams.ndof);

for c=1:rodParams.ne
    node_1 = getVertex(rodParams.x, sElement(c).nodeIndex(1));
    node_2 = getVertex(rodParams.x, sElement(c).nodeIndex(2));

    [dF, dJ] = stretchingForce(node_1, node_2, sElement(c).refLen, sElement(c).EA_local);
    
    index = sElement(c).globalIndex;
    
    Fs(index) = Fs(index) - dF;
    Js(index,index) = Js(index,index) - dJ;
end

end
