% -- Set up some variables which define the optimisation

% =================================================================
%    IT'S MOSTLY UP TO YOU TO MAKE SURE THESE THINGS ARE COMPATIBLE
% =================================================================

nSide = 4;                   % -- the number of vertices of your polygon
symOrb = [1 4 1];            % -- the number of points in each of the three symmetry orbits
nTrial = 100;                % -- the number of trial bases to use
realBasisType = 'Maximal2D'; % -- the basis to be used over the solution points
testBasisType = 'Maximal2D'; % -- the test basis to which ours is compared
wT = ones(nTrial,1);         % -- the weighting of the terms in the test basis
p = 2;                       % -- the p-norm magnitude
nPop = 5000;                 % -- the GA optimisation population size
nGen = 100;                  % -- the GA optimisation maximum number of generations

% =========================================================================
%    YOU MAY WELL FIND THAT SOME SYMMETRY ORBITS AND BASES ARE INCOMPATIBLE
% =========================================================================

% -- Add the path to the basis functions subroutine to keep it all in one place
addpath('../basisFunctions/');

% -- Add the path to the routines for finding quadratures on polygons to keep it all in one place
addpath('../integrationWeights/');

% -- Set up the polygon and the symmetry matrix
[refPolygon, nSoln, nOutr, xOutr, symMatx] = getPolygonSymmetry(nSide, symOrb);
fprintf('Symmetry generating functions calculated.\n')

% -- Set up the integration points and weights over the polygon
[xIntg, wIntg] = getIntegrationPoints(refPolygon, 20);
fprintf('Integration points and weights over polygon found.\n')

% -- Set up the actual and trial basis functions over the polynomial, and orthonormalise
realBasisN = getBasisFunctions(realBasisType, nSoln,  refPolygon);
testBasisN = getBasisFunctions(testBasisType, nTrial, refPolygon);
fprintf('Basis functions set and orthonormalised.\n')

% -- Set up the bounds of the optimisation
[lb, ub] = getOptimisationBounds(symOrb);
fprintf('Bounds of optimisation set up.\n')

% -- Set up the function(s) to be optimised
getError = @(oVec) getPointSetError(getPointSetFromVector(symMatx, symOrb, oVec), realBasisN, testBasisN, xIntg, wIntg, wT, p);
fprintf('Optimisation error function set.\n')

% -- Set up the optimisation settings (genetic algorithm)
options = optimoptions('ga','Display','iter', ...
                            'MaxStallGenerations',floor(nGen/3), ...
                            'MaxGenerations',nGen, ...
                            'CrossoverFcn',{@crossoverintermediate}, ...
                            'OutputFcn',{@(options,state,flag) gaoutfun(options, state, flag, nOutr, symOrb, symMatx)}, ...
                            'PopulationSize',nPop, ...
                            'PlotFcn',{@gaplotbestf, ...
                                       @gaplotscorediversity}, ...
                            'UseParallel',true);
fprintf('Optimisation settings defined.\n')

% -- Perform the optimisation (genetic algorithm)
[sVecOpt, fVal] = ga(getError, length(lb), [], [], [], [], lb, ub, [], options);

% ----------------------------------------------------
% OPTIMISATION VISUALISATION AND OUTPUT FUNCTIONS
% ----------------------------------------------------

% -- Output function for single target Genetic Optimisation

function [state,options,optchanged] = gaoutfun(options, state, flag, nSide, sOrb, SymMat)
optchanged = false;
switch flag
    case 'iter'

        % -- Update the history every 50 generations.
        if rem(state.Generation,50) == 0

            % -- Find the best objective function, and stop if it is low.
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);

            % -- Print the best performing set to screen for "safekeeping"
            xT = getPointSetFromVector(SymMat, sOrb, bestx)

        end
end

end
