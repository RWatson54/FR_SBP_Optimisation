function [SBPError, CSVError] = getSBP_Numeric(xOutr, xShap, xSoln, xFlux, nFluF, shapBasis, compBasis, Integrate)

% -- This code tests the SBP property of the matrix kernels **numerically** for 2D FR

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

%       xOutr -- the coordinates of the outside points, useful for normals and integration area (anticlockwise numbering, must be convex)
%       xShap -- the coordinates of the points which will be used to geometrically fit the shape (so, might have 9 points for a P2 quad)
%       xSoln -- the coordinates of the solution points within the shape
%       xFlux -- the coordinates of the flux points, compatible with the order defined by nFluF
%       nFluF -- the number of flux points on each face (e.g., [5 1 5 1]). Must be compatible with xFlux and xOutr
%       shapBasis -- a function handle to return the evaluated basis to be used for geometrically fitting the shape at coordinates, as a function of x and d [= @(x,d) x(1) + x(2);]
%       compBasis -- a function handle to return the evaluated basis to be used for determining the solution computationally, as a function of x and d [= @(x,d) x(1) + x(2);]

% -- And the outputs are:
%       SBPError -- the Summation-by-Parts error of the scheme
%       CSVError -- the conservation error of the scheme

% -- Get the number of dimensions from the outer points
nDim = size(xOutr,2);

% -- Set up the computational domain, useful for plotting easily
compPoly = polyshape(xOutr);
% fprintf('Computational domain set\n')

% -- The integration points
nlInt = 200; % The number of line integration points per face
if ~exist('Integrate','var') || isempty(Integrate)
    % If the integration technique isn't specified, set it to zero
    nInt = 100;
    [xIntegrate, wIntegrate] = getIntegrationPoints(compPoly,nInt);
else
    if isa(Integrate, 'double')
        % If the variable Integrate contains a number, then use this as the integration weight per triangle 
        nInt = Integrate;
        [xIntegrate, wIntegrate] = getIntegrationPoints(compPoly,nInt);
    elseif (isa(Integrate, 'char') || isa(Integrate, 'string'))
        % If the input is a string or character vector, then treat this as a filename stored in "../integrationWeights/explicit/<<filename>>"
        dataIn = readmatrix(strjoin(["../integrationWeights/explicit/", Integrate],''));
        nInt = size(dataIn,1); xIntegrate = dataIn(:,1:2); wIntegrate = dataIn(:,3);
    else
        % Otherwise, no idea what to do
        error(' Unsure how to perform cubature over the element, terminating ')
    end

end

% -- Extract the surface integration points and corresponding normals from xOutr
nFace = size(nFluF,2);
for iFace = 1:nFace

    % -- Work out which of the outer points bound this face
    lP = mod(iFace  -1,nFace)+1;
    rP = mod(iFace+1-1,nFace)+1;

    % -- Get the direction of the normal vector
    n(1) = (xOutr(rP,2) - xOutr(lP,2));
    n(2) = (xOutr(lP,1) - xOutr(rP,1));

    % -- Get the integration weights along the edge
    [xFInt((iFace-1)*nlInt+1:iFace*nlInt,:), wFInt((iFace-1)*nlInt+1:iFace*nlInt,:)] = lgwt(nlInt, xOutr(lP,:), xOutr(rP,:));

    % -- And add in the normalised normals
    xFNrm((iFace-1)*nlInt+1:iFace*nlInt,1) = n(1) / norm(n);
    xFNrm((iFace-1)*nlInt+1:iFace*nlInt,2) = n(2) / norm(n);

end
% fprintf('Integration points, weights, and normals successfully set\n')

% -- And then compute the orthonormalised bases
sBasisSet = shapBasis(xIntegrate,0,xShap);
gramShapO = ((sBasisSet .* wIntegrate') * sBasisSet');
cBasisSet = compBasis(xIntegrate,0,xSoln);
gramCompO = ((cBasisSet .* wIntegrate') * cBasisSet');
% fprintf('Gramian matrices calculated\n')

% -- By using the Cholesky decomposition of the Gramian as a quick MGS
basisShapN = @(x,d,xC) inv(chol(gramShapO))' * shapBasis(x,d,xC);
basisCompN = @(x,d,xC) inv(chol(gramCompO))' * compBasis(x,d,xC);
% fprintf('Basis functions orthonormalised over domain\n')

% -- Plot the domain and the flux and solution points, for a check
%plot(compPoly); hold on; plot(xSoln(:,1), xSoln(:,2),'x'), plot(xFlux(:,1), xFlux(:,2), 'o')
% fprintf('Flux and solution point locations set\n')

% -- Start doing the work by building some Alternants
AlternantCSO = basisCompN(xSoln,0,xSoln);
AlternantCQO = basisCompN(xIntegrate,0,xSoln);
AlternantCFO = basisCompN(xFlux,0,xSoln);
AlternantCLO = basisCompN(xFInt,0,xSoln);
for iDim = 1:nDim
    AlternantCSD{iDim} = basisCompN(xSoln,iDim,xSoln);
end
% fprintf('Alternant matrices formed on the computational basis\n')

% -- Get the mass matrix for the quadrature at the solution points
% -- diag(Msq) .* M === diag(wIntegrate) * M, so could leave as a vector
Msq = diag(wIntegrate);
% fprintf('Quadrature mass matrix formed\n')

% -- Set up the "augmented mass matrix", for use with the gradients
% -- diag(Msq_d) .* M === kron(eye(2), diag(wIntegrate)) * M, so could leave as a vector
Msq_d = kron(eye(2), diag(wIntegrate));
% fprintf('Augmented quadrature mass matrix formed\n')

% -- Get the solution-to-quadrature projection matrix
Lsq = (AlternantCSO \ AlternantCQO)';
% fprintf('Solution-to-quadrature projection matrix formed\n')

% -- Set up the "augmented solution-to-quadrature projection matrix
Lsq_d = kron(eye(nDim), Lsq);
% fprintf('Augmented Solution-to-quadrature projection matrix formed\n')

% -- Build the differentiation matrices
for iDim = 1:nDim
    M2t{:,iDim} = cleanZeros(double((AlternantCSO) \ (AlternantCSD{iDim}))');
end
M2 = cell2mat(M2t);
M5 = cell2mat(M2t');

% -- Matrix D is equivalent to the divergence operator, M2
D = M2;
% fprintf('Solution point flux divergence matrix found\n')

% -- Matrix G is equivalent to the differentiation operator, M5
G = M5;
% fprintf('Solution point primitive differentiation matrix found\n')

% -- Assemble the LHS of the SBP property
SBP_lhs = (Lsq'*Msq*Lsq*D) + (G'*Lsq_d'*Msq_d*Lsq_d);

% -- Matrix LSF is equivalent to the solution to flux point projection operation, M1
M1 = cleanZeros(double(((AlternantCSO) \ (AlternantCFO))'));
Lsf = M1;
% fprintf('Solution point to flux point primitive projection matrix found\n')

% -- Set up the "augmented" solution-to-flux projection matrix
Lsf_d = kron(eye(nDim), Lsf);
% fprintf('Augmented solution point to flux point projection matrix formed\n')

% -- Get the boundary mass matrix for the quadrature at the flux points
% -- diag(Mfq) .* M === diag(wFInt) * M, so could leave as a vector
Mfq = diag(wFInt);
% fprintf('Boundary quadrature mass matrix formed\n')

% -- Matrix LFR isn't equivalent to any normal FR matrix!
for iFace = 1:nFace

    % -- It's pretty frustrating to try and edit a line this long, so break it down
    ir1 = (iFace-1)*nlInt+1:iFace*nlInt;
    ir2 = sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace));

    Lfr(ir1,ir2) = cleanZeros(double((AlternantCFO(:,ir2) \ AlternantCLO(:,ir1))'));

end
%Lfr = cleanZeros(double(((AlternantCFO) \ (AlternantCLO))'));
% fprintf('Flux point to boundary quadrature point projection matrix found\n')

% -- Set up the "augmented" flux to boundary quadrature point projection matrix
Lfr_d = kron(eye(nDim), Lfr);
% fprintf('Augmented flux to boundary quadrature point projection matrix formed\n')

% -- Put the normals into the preferred form for matrix multiplication
N = zeros(size(xFNrm,1), nDim*size(xFNrm,1));
for iDim = 1:nDim
    N(:,(iDim-1)*size(xFNrm,1)+1:iDim*size(xFNrm,1)) = diag(xFNrm(:,iDim));
end
% fprintf('Matrix of normals formed\n')

% -- Assemble the right hand side
SBP_rhs = Lsf'*Lfr'*Mfq*N*Lfr_d*Lsf_d;

% -- And determine the SBP Property
SBPError = norm(SBP_lhs - SBP_rhs);
% fprintf('SBP Property determined\n')

% ------------------------------------------------------------------------- %
% ----------------------- SBP PROPERTY DETERMINED ------------------------- %
% ------------------------------------------------------------------------- %

% -- The C matrix is just the correction function matrix, M4
for iFace = 1:nFace

    % -- Work out which of the outer points bound this face
    lP = mod(iFace  -1,nFace)+1;
    rP = mod(iFace+1-1,nFace)+1;

    % -- Get the integration weights along the edge
    [x1DInt, w1DInt] = lgwt(nlInt, xOutr(lP,:), xOutr(rP,:));

    % -- Use the MP pseudoinverse to solve the "overconstrained" system for each face
    lT = pinv(basisCompN(xFlux(sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace)),:),0,xSoln))*basisCompN(x1DInt,0,xSoln);

    % -- Also evaluate the orthonormal bases themselves along the line
    oT = basisCompN(x1DInt,0,xSoln);

    % -- Numerically integrate along the line to give the coefficients
    sTC(:,sum(nFluF(1:iFace-1))+1:sum(nFluF(1:iFace))) = ((oT .* w1DInt') * lT');

end
M4 = cleanZeros(double(AlternantCSO' * sTC));
C = M4;
% fprintf('Correction function matrix found\n')

% -- Assemble the left hand side
CSV_lhs = ones(size(xSoln,1),1)'*Lsq'*Msq*Lsq*C;

% -- Assemble the right hand side
CSV_rhs = ones(size(xSoln,1),1)'*Lsf'*Lfr'*Mfq*Lfr;

% -- And calculate the norm 
CSVError = norm(CSV_lhs - CSV_rhs);

end

% -- HELPER FUNCTIONS
% -- CleanZeros: tidy up very small numbers from matrices
function M = cleanZeros(M)

    % -- This function just cleans up any very small values

    % -- Set the tolerance
    tol = 0.0000000001;

    % Zero the values which fall below this
    M(abs(M)<tol) = 0; 

end
% -- Interleave: thread matrices together row by row / column by column
function output = interleave(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function interleaves the rows or columns of the specified vectors or
% matrices A, B, C, .... The last argument indicates whether to interleave
% the rows or columns by passing 'row' or 'col', respectively. If
% interleaving the rows, the number of columns of A and B must be the same;
% if interleaving the columns, the number of rows must be the same. 
%
% The order of input vectors/matrices indicate the order of interleaving
% (e.g. first row of A will be the first row, first row of B will be the
% second row of the interleaved matrix).
% 
% Extra rows or columns will be appended to the end of the final vector or
% matrix.
%
% If all inputs are vectors, their elements will be interleaved without
% the vectors needing to be all rows or all columns if the last argument
% ('row' or 'col') is not specified. If the orientation of the interleaving
% is specified, then the program will follow that orientation.
% 
% Examples:
% 1) Interleaving rows of matrices
% A = [1 2           B = [5 6
%      3 4]               7 8]     
% C = interleave(A, B, 'row')
% C = [1 2
%      5 6
%      3 4
%      7 8]
%
% 2) Interleaving columns of matrices
% C = interleave(A, B, 'col')
% C = [1 5 2 6
%      3 7 4 8]
%
% 3) Interleaving vectors
% A = [1 2 3 4]      B = [5 6 7 8 9]'
% C = interleave(A, B)
% C = [1 5 2 6 3 7 4 8 9]
%
% 4) Interleaving >2 matrices
% A = [1 2           B = [5 6         C = [9 10      D = [13 14
%      3 4]               7 8]            11 12]          15 16]
% E = interleave(A, B, C, D, 'col')
% E = [1 5  9 13 2 6 10 14
%      3 7 11 15 4 8 12 16]
% 
% 5) Interleaving columns of 2 matrices with unequal columns
% A = [1 2           B = [5 6  7  8
%      3 4]               9 10 11 12] 
% C = interleave(A, B, 'col')
% C = [1 5 2 6  7  8
%      3 9 4 10 11 12]
%     
% 6) Interleaving >2 vectors of unequal lengths
% A = [1 2 3 4]    B = [5 6 7]    C = [8 9 10 11 12 13]
% D = interleave(A, B, C, D)
% D = [1 5 8 2 6 9 3 7 10 4 11 12 13]
% 
% Written by Steven Shimizu
% March 4, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if interleaving orientation is specified
if ischar(varargin{end}) %check if last argument is a char array
    orientation = varargin{end};
    
    %Then check if interleaving orientation is anything but 'row' or 'col'
    if ~(isequal(orientation,'row') || isequal(orientation,'col'))
        %if not equal to 'row' or 'col', throw an error
        error('Last argument should specify ''row'' or ''col'' unless the inputs only consist of vectors.');
    else
        %The orientation is specified, then there are n-1 input numeric
        %arrays (single numbers, vectors or 2d matrices).
        numInputs = nargin - 1;
        
        for i = 1:1:numInputs
            if ~isnumeric(varargin{i})
                error('Inputs must be numeric arrays.');
            end
        end
    end
    
%check if last input argument is a numeric array, then check that all
%of the input arguments are vectors only, otherwise throw an error.
elseif isnumeric(varargin{end})
    orientation = 'none';
    numInputs = nargin;
    for i = 1:1:numInputs
        if ~isvector(varargin{i})
            error('The interleaving orientation ''row'' or ''col'' was not specified, which is only valid if all inputs are vectors.');
        end
    end
    
end
%Now that all inputs are known to be numeric arrays and interleaving
%orientation is specified, check that they are the right size
sizeArray = zeros(numInputs,2);
lengthList = zeros(numInputs,1);
for i = 1:1:numInputs
    sizeArray(i,:) = size(varargin{i});
    if isequal(orientation,'none')
        lengthList(i) = length(varargin{i});
    end   
end
if isequal(orientation,'row') && ~(sum(sizeArray(:,2)==sizeArray(1,2)) == numInputs)
    %check if number of columns are the same
    error('The number of columns of input matrices must be the same for interleaving the rows');
elseif isequal(orientation,'col') && ~(sum(sizeArray(:,1)==sizeArray(1,1)) == numInputs)
    error('The number of rows of input matrices must be the same for interleaving the columns');
end
output = [];        
%If only a single numeric array is passed, then return the same array
if numInputs == 1
    output = varargin{1};
elseif isequal(orientation,'none')
    %interleave vectors
    output = zeros(sum(lengthList),1);
    maxLength = max(lengthList);
    i_output = 1;
    for i = 1:1:maxLength
        for n = 1:1:numInputs
            if i <= lengthList(n)
                output(i_output) = varargin{n}(i);
                i_output = i_output + 1;
            end
        end
    end
elseif isequal(orientation,'row')
    %interleave by rows
    output = zeros(sum(sizeArray(:,1)),sizeArray(1,2));
    maxRows = max(sizeArray(:,1));
    i_output = 1;
    for i = 1:1:maxRows
        for n = 1:1:numInputs
            if i <= sizeArray(n,1)
                output(i_output,:) = varargin{n}(i,:);
                i_output = i_output + 1;
            end
        end
    end
elseif isequal(orientation,'col')
    %interleave by columns
    output = zeros(sizeArray(1,1),sum(sizeArray(:,2)));
    maxCols = max(sizeArray(:,2));
    i_output = 1;
    for i = 1:1:maxCols
        for n = 1:1:numInputs
            if i <= sizeArray(n,2)
                output(:,i_output) = varargin{n}(:,i);
                i_output = i_output + 1;
            end
        end
    end
end
end