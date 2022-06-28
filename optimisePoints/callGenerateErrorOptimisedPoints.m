% -- Set up some variables which define the optimisation

% -- Tidy up - you know there are going to be problems if you don't
close all
clear all

% =================================================================
%    IT'S MOSTLY UP TO YOU TO MAKE SURE THESE THINGS ARE COMPATIBLE
% =================================================================

nSide = 4;                         % -- the number of vertices of your polygon
symOrb = [1 4 1];                  % -- the number of points in each of the three symmetry orbits
nFace = 5;                          % -- the number of points per side along each face (Gauss-Legendre distribution)
nTrial = 144;                      % -- the number of trial bases to use
compBasisType.Type = 'Maximal2D'; % -- the basis to be used over the solution points
compBasisType.Eps  = 1e-9;
compBasisType.xC   = 'Dummy';
testBasisType.Type = 'Maximal2D';  % -- the test basis to which ours is compared - strongly recommend polynomial here, since avoids the need to lay out centres over the domain
testBasisType.Eps  = 'Dummy';
testBasisType.xC   = 'Dummy';
wT = ones(nTrial,1);               % -- the weighting of the terms in the test basis
p = 2;                             % -- the p-norm magnitude
nPop = 500;                       % -- the GA optimisation population size
nGen = 25;                        % -- the GA optimisation maximum number of generations

% =================================================================
%    IT'S MOSTLY UP TO YOU TO MAKE SURE THESE THINGS ARE COMPATIBLE
% =================================================================

% -- Call the optimisation
[sVecOpt, fVal] = generateErrorOptimisedPoints(nSide, symOrb, nFace, nTrial, compBasisType, testBasisType, wT, p, nPop, nGen);
