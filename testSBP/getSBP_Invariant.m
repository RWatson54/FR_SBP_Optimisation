function [Msq, Msq_d, Mfq, N] = getSBP_Invariant(wIntI, wIntF, nIntF)
%GETSBP_EXTERNAL Summary of this function goes here

% -- Get the dimensions from the coordinates of the normals
nDim = size(nIntF,2);

% -- First, get the mass matrix over the solution points
Msq = diag(wIntI);

% -- Set up the "augmented mass matrix", for use with the gradients
Msq_d = kron(eye(2), Msq);

% -- Get the boundary mass matrix for the quadrature at the flux points
Mfq = diag(wIntF);

% -- Put the normals into the preferred form for matrix multiplication
N = zeros(size(nIntF,1), nDim*size(nIntF,1));
for iDim = 1:nDim
    N(:,(iDim-1)*size(nIntF,1)+1:iDim*size(nIntF,1)) = diag(nIntF(:,iDim));
end

end

