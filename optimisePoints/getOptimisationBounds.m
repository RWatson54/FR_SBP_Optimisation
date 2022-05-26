function [lb, ub] = getOptimisationBounds(symOrb)

%GETOPTMISATIONBOUNDS generates the bounds on each value of the symmetry orbits needed by the optimiser

% -------------------------------------------------------------

% -- The inputs are:

%       symOrb -- the enumerated symmetry orbits, which effects the bounds

% -- The outputs are:

%       lb -- the lower bounds of the optimisation
%       ub -- the upper bounds of the optimisation

% -------------------------------------------------------------

% -- Zero the counter
iX = 0;

% -- No bounds on the initial value, because it isn't really there...
for iOrb = 1:symOrb(1)
end

% -- For the single-point ones, these can vary between -1 and 1
for iOrb = 1:symOrb(2)
    iX = iX + 1;
    lb(iX) = -1; ub(iX) = 1;
end

% -- The two pointers vary between 0 and 1
for iOrb = 1:symOrb(3)
    iX = iX + 2;
    lb(iX-1) = 0; ub(iX-1) = 1;
    lb(iX  ) = 0; ub(iX  ) = 1;
end

end
