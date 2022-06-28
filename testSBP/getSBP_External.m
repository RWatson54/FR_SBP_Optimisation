function [SBPError] = getSBP_External(xOutr, xSoln, xFlux, xIntI, xIntF, nFluF, compBasis, Msq, Msq_d, Mfq, N)
%GETSBP_EXTERAL Summary of this function goes here

% -- This code tests the SBP property of the matrix kernels **numerically** for 2D FR, but takes a lot of the work elsewhere

% ============================================================================== %
%
%    I felt like I should probably warn you here - this is my first experience with 
%    SBP, and I wasn't completely sure what I was doing. It might be worth checking
%    this over quickly to make sure it makes sense before relying heavily on it...
%
%       -- Rob
%
% ============================================================================== %

% -- The inputs are:

%       xOutr -- the coordinates of the outer points of the shape
%       xSoln -- the coordinates of the solution points within the shape
%       xFlux -- the coordinates of the flux points, compatible with the order defined by nFluF
%       xIntI -- the integration points over the interior of the element
%       xIntF -- the integration points over the faces of the element
%       Msq -- the mass matrix over the interior of the element
%       Msq_d -- the augmented mass matrix
%       Mfq -- the mass matrix over the faces of the element
%       N -- the matrix of normals
%       compBasis -- a function handle to return the evaluated basis to be used for determining the solution computationally, as a function of x and d [= @(x,d) x(1) + x(2);]

% -- And the outputs are:
%       SBPError -- the Summation-by-Parts error of the scheme

% ============================================================================== %

% -- Get the number of dimensions
nDim = size(xSoln,2);

% -- Start doing the work by building some Alternants
AlternantCSO = compBasis(xSoln,0,xSoln);
AlternantCQO = compBasis(xIntI,0,xSoln);
AlternantCFO = compBasis(xFlux,0,xSoln);
AlternantCLO = compBasis(xIntF,0,xSoln);
for iDim = 1:nDim
    AlternantCSD{iDim} = compBasis(xSoln,iDim,xSoln);
end

% -- Get the solution-to-quadrature projection matrix
Lsq = (AlternantCSO \ AlternantCQO)';

% -- Set up the "augmented solution-to-quadrature projection matrix
Lsq_d = kron(eye(nDim), Lsq);

% -- Build the differentiation matrices
for iDim = 1:nDim
    M2t{:,iDim} = (AlternantCSO \ AlternantCSD{iDim})';
end
D = cell2mat(M2t);
G = cell2mat(M2t');

% -- Assemble the LHS of the SBP property
SBP_lhs = (Lsq'*Msq*Lsq*D) + (G'*Lsq_d'*Msq_d*Lsq_d);

% -- Matrix LSF is equivalent to the solution to flux point projection operation, M1
Lsf = (AlternantCSO \ AlternantCFO)';

% -- Set up the "augmented" solution-to-flux projection matrix
Lsf_d = kron(eye(nDim), Lsf);

% -- Matrix LFR isn't equivalent to any normal FR matrix!
nlInt = size(xIntF,1) / size(xOutr,1);
for iFace = 1:size(xOutr,1)

    % -- It's pretty frustrating to try and edit a line this long, so break it down
    ir1 = (iFace-1)*nlInt+1:iFace*nlInt;
    ir2 = sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace));

    Lfr(ir1,ir2) = (AlternantCFO(:,ir2) \ AlternantCLO(:,ir1))';

end

% -- Set up the "augmented" flux to boundary quadrature point projection matrix
Lfr_d = kron(eye(nDim), Lfr);

% -- Assemble the right hand side
SBP_rhs = Lsf'*Lfr'*Mfq*N*Lfr_d*Lsf_d;

% -- And determine the SBP Property
SBPError = norm(SBP_lhs - SBP_rhs);

end

