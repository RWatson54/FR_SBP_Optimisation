% -- This is a driving script for generateFRMatrixKernels
% -- Because it's numerical, it relies on the getBasisFunctions script in ../basisFunctions

% -- Tidy up from a previous run
close all
clear all

% -- Set the number of points
nOutr = 4; nFace = 5;

% -- Add the path to the basis functions subroutine to keep it all in one place
addpath('../basisFunctions/');

% -- Add the path to the integration weights subroutines to keep them all in one place
addpath('../integrationWeights/');

% -- Set the coordinates of the outer points which define the shape of the element
compPoly = nsidedpoly(nOutr,'Center',[0 0],'SideLength',2);
xOutr = flipud(compPoly.Vertices);

% -- Set the coordinates of the "shape" points which define the control points for positioning and scaling the element
xShap = xOutr;

% -- Set the solution point coordinates required, could put in directly, could Blue Peter read from text file
xSoln = readmatrix('../optimisePoints/pointSets/solnSet_d2_p4_n25_001.txt');
xSoln = xSoln + 2*0.00001*(randn(size(xSoln,1), size(xSoln,2)) - 0.5);

% -- Set up the distribution of points around the element
nFluF = repmat(nFace, 1, nOutr);

% -- Set up the computational flux points
% Flux Points
% Set up the flux point "orbits" - one orbit for each side
for iF = 1:nOutr
    fluxOrbit(iF).Side = lgwt(nFace,-1,1)'; % -- Gauss points
end
xFlux = [];
for iF = 1:size(fluxOrbit,2)
    
    % Get the vector pointing along the edge
    dx = xOutr(mod(iF,size(xOutr,1))+1,:) - xOutr(iF,:);
    
    for iP = 1:size(fluxOrbit(iF).Side,2)
        
        % Get the new set of flux points
        xT = [(0.5 * fluxOrbit(iF).Side(iP) + 0.5) * dx + xOutr(iF,:)];
        xFlux = [xFlux; xT];

    end
end

% -- Set the basis type for the shape and computational basis
compBasisType.Type = 'GaussianGA';
compBasisType.xC = xSoln;
compBasisType.Eps = 1e-10;

% -- Set the function handles for the bases
compBasis = getBasisFunctions(compBasisType,size(xSoln,1));

% -- Call the generating code, looping over various quadrature schemes a la Trojak

quadraturesList = dir('../integrationWeights/explicit/poly_m6*.txt');
quadFilenames = string({quadraturesList.name});

% -- Sort things out
[xIntI, wIntI] = getIntegrationPoints(compPoly, 80);
[xIntF, wIntF, nIntF] = getFaceIntegration(xOutr, 200);
[Msq, Msq_d, Mfq, N] = getSBP_Invariant(wIntI, wIntF, nIntF);

SBP = getSBP_External(xOutr, xSoln, xFlux, xIntI, xIntF, nFluF, compBasis, Msq, Msq_d, Mfq, N)
