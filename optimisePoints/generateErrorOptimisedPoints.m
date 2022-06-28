function [xOut, fVal] = generateErrorOptimisedPoints(nSide, symOrb, nFace, nTrial, compBasisType, testBasisType, wT, p, nPop, nGen)

% =========================================================================
%    YOU MAY WELL FIND THAT SOME SYMMETRY ORBITS AND BASES ARE INCOMPATIBLE
% =========================================================================

% -- Add the path to the basis functions subroutine to keep it all in one place
addpath('../basisFunctions/');

% -- Add the path to the routines for finding quadratures on polygons to keep it all in one place
addpath('../integrationWeights/');

% -- Add the path to the SBP tester to keep things in one place
addpath('../testSBP/');

% -- Set up the polygon and the symmetry matrix
[refPolygon, nSoln, nOutr, xOutr, xFlux, nFluF, symMatx] = getPolygonSymmetry(nSide, symOrb, nFace);
fprintf('Symmetry generating functions calculated.\n')

% -- Set up the integration points and weights over the polygon interior
[xIntI, wIntI] = getIntegrationPoints(refPolygon, 40);
fprintf('Integration points and weights over polygon interior found.\n')

% -- Set up the integration points and weights over the polygon faces
[xIntF, wIntF, nIntF] = getFaceIntegration(xOutr, 50);
fprintf('Integration points and weights over polygon faces found.\n')

% -- Set up the various invariant matrices needed by the SBP test
[Msq, Msq_d, Mfq, N] = getSBP_Invariant(wIntI, wIntF, nIntF);

% -- Set up the actual and trial basis functions over the polynomial, and orthonormalise
compBasisN = getBasisFunctions(compBasisType, nSoln);
testBasisN = getBasisFunctions(testBasisType, nTrial, refPolygon);
fprintf('Basis functions set and test function orthonormalised.\n')

% -- Set up the bounds of the optimisation
[lb, ub] = getOptimisationBounds(symOrb);
fprintf('Bounds of optimisation set up.\n')

% -- Set up the function(s) to be optimised
getError = @(oVec) getPointSetError(xOutr, xFlux, getPointSetFromVector(symMatx, symOrb, oVec), compBasisN, testBasisN, xIntI, wIntI, xIntF, nFluF, Msq, Msq_d, Mfq, N, wT, p);
fprintf('Optimisation error function set.\n')

% -- Set up the optimisation settings (genetic algorithm)
options = optimoptions('ga','Display','iter', ...
                            'MaxStallGenerations',min(floor(nGen/5), 50), ...
                            'MaxGenerations',nGen, ...
                            'CrossoverFcn',{@crossoverintermediate}, ...
                            'OutputFcn',{@(options,state,flag) gaoutfun(options, state, flag, nOutr, symOrb, symMatx)}, ...
                            'PopulationSize',nPop, ...
                            'UseParallel',true);
                            %'PlotFcn',{@gaplotbestf, ...
                            %           @gaplotscorediversity}, );
fprintf('Optimisation settings defined.\n')

% -- Perform the optimisation (genetic algorithm)
[sVecOpt, fVal] = ga(getError, length(lb), [], [], [], [], lb, ub, [], options);

xOut = getPointSetFromVector(symMatx, symOrb, sVecOpt);

end

% ----------------------------------------------------
% OPTIMISATION VISUALISATION AND OUTPUT FUNCTIONS
% ----------------------------------------------------

% -- Output function for single target Genetic Optimisation

function [state,options,optchanged] = gaoutfun(options, state, flag, nSide, sOrb, SymMat)
optchanged = false;
switch flag
    case 'iter'

        % -- Update the history every 50 generations.
        if rem(state.Generation,25) == 0

            % -- Find the best objective function, and stop if it is low.
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);

            % -- Print the best performing set to screen for "safekeeping"
            xT = getPointSetFromVector(SymMat, sOrb, bestx)

        end
end

end
