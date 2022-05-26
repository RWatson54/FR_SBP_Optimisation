function [basis] = getBasisFunctions(basisType, nBasis, compPoly)

%GETBASISFUNCTIONS generates a function handle which evaluates the requested nBasis bases as a column vector at each of the input coordinates

% -------------------------------------------------------------

% -- The inputs are:

%       basisType -- The type of the basis that is being requested
%       nBasis -- how many bases of this type are being asked for
%       compPoly (optional) -- the polygon that we're using here, for the orthonormalisation

% -- The outputs are:

%       basis -- the function handle of the orthonormalised basis

% -------------------------------------------------------------

% -- Test to see if the basis number requested is reasonable
switch basisType

    case 'Maximal2D'

        % -- For the 2D maximal order basis, the requested number of bases should be a square number
        isSquare = mod(sqrt(nBasis),1)==0;
        if ~isSquare

            error('Incompatible number of bases requested for a Maximal2D basis')

        end

    otherwise

        error('No corresponding trial basis functions set. Goodbye.\n')

end

% -- Set up the bases for the requested type
switch basisType

    case 'Maximal2D'

        % -- For the 2D maximal order basis, the requested number of bases should be a square number
        basisO = @(x,d) maximalOrder2D(x,d,nBasis);

    otherwise

        error('No corresponding trial basis functions set. Goodbye.\n')

end

% -- Test to see if the polygon shape has been passed in for orthonormalisation
if exist('compPoly','var')
    % -- Generate the integration weights
    [xInt, wInt] = getIntegrationPoints(compPoly, 20);
    % -- Integrate the basis over the domain to get the Gramian
    Gm = (basisO(xInt,0) .* wInt') * basisO(xInt,0)';
    % -- And calculate the matrix for orthonormalisation with the inverse of the Cholesky
    oMat = inv(chol(Gm))';
else
    % -- Otherwise, set the orthonormalisation matrix to an identity
    oMat = eye(nBasis);
end

% -- Othonormalise the bases if required
basis = @(x,d) oMat * basisO(x,d);

end

% ---------------------------------------------
% ---------------------------------------------
%  Basis generating functions
% ---------------------------------------------
% ---------------------------------------------

% Feel free to add your own basis functions in here in the style of maximalOrder2D!

% -- Function for two dimensional maximal order basis
function [b] = maximalOrder2D(x,d,n)

    % -------------------------------------------------------------

    % -- The inputs are:

    %       x -- the coordinates to evaluate the basis at [nCoord -x- nDim]
    %       d -- the direction in which to differentiate, if any
    %       n -- the number of bases to compute

    % -- The outputs are:

    %       b -- the basis evaluated at the coordinate points

    % -------------------------------------------------------------

    % -- Get the powers involved
    nOMaxOrder = round(sqrt(n));
    pow1D = 0:1:nOMaxOrder-1;
    xPow = reshape(repmat(pow1D',1,nOMaxOrder),[],1);
    yPow = reshape(repmat(pow1D ,nOMaxOrder,1),[],1);

    % -- Compute the basis
    if d == 0
        b = (x(:,1)'.^(xPow) .* x(:,2)'.^(yPow));
    elseif d == 1   
        b = (max(0,xPow) .* x(:,1)'.^(max(0,xPow-1)) .* x(:,2)'.^(yPow));   
    elseif d == 2
        b = (x(:,1)'.^(xPow) .* max(0,yPow) .* x(:,2)'.^(max(0,yPow-1)));
    end

end

