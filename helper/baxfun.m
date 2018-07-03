function c = baxfun(f, a, b, shiftA, shiftB)
% BAXFUN  Binary array expansion function
%   C = BAXFUN(FUNC, A, B) applies an element-by-element binary operation
%   to arrays A and B, with singleton expansion (SX) enabled. It is
%   equivalent to C = BSXFUN(FUNC, A, B). FUNC is a function handle. 
% 
%   C = BAXFUN(FUNC, A, B, SHIFTA, SHIFTB) applies an element-by-element
%   binary operation to arrays A and B, with array expansion (AX) enabled.
%   FUNC is a function handle. SHIFTA and SHIFTB specify the number of
%   dimensions by which arrays A and B are optionally shifted to the right
%   (positive values) or left (negative values). Thus, 
%   C = BAXFUN(FUNC, A, B, SHIFTA, SHIFTB) is equivalent to 
%   C = BAXFUN(FUNC, SHIFTDIM(A,-SHIFTA), SHIFTDIM(B,-SHIFTB)).
%   C = BAXFUN(FUNC, A, B, SHIFTA) can be used if B is not to be shifted.
% 
%   FUNC can either be an handle for an M-function, or one of the following
%   built-in-function handles:
%
%       @plus     Plus                     @eq    Equal
%       @minus    Minus                    @ne    Not equal
%       @times    Array multiply           @lt    Less than
%       @rdivide  Right array divide       @le    Less than or equal
%       @ldivide  Left array divide        @gt    Greater than
%       @power    Array power              @ge    Greater than or equal
%       @max      Binary maximum           @and   Logical AND
%       @min      Binary minimum           @or    Logical OR
%       @rem      Division remainder       @xor   Logical EXCLUSIVE OR
%       @mod      Division modulus
%       @atan2    4-quadrant arc-tangent
%       @hypot    SQRT of sum of squares
%
%   If an M-function handle is specified, the M-function must be able to
%   accept as input either two column vectors of the same size, or one
%   column vector and one scalar, and return as output a column vector of
%   the same size as the input(s). 
%
%   Array expansion (AX) is a powerful generalization of the concept of
%   scalar expansion. Scalar expansion is the virtual replication or
%   annihilation of a scalar which allows you to combine it, element by
%   element, with an array X of any size (e.g. X+10, X*10, or []-10).
%   Similarly, in BAXFUN, the purpose of AX is to virtually match the sizes
%   of A and B, before FUNC is applied. Indeed, A and B can be scalars,
%   vectors, matrices, or multi-dimensional arrays, and may have different
%   sizes. Dimension matching is achieved by means of a dimension shift
%   followed by a singleton expansion:
% 
%   1) DIMENSION SHIFT
%      A dimension shift (see SHIFTIM) is applied to A or B if a non-zero
%      value is specified for SHIFTA or SHIFTB. Notice that, in BAXFUN, a
%      shift to the left is hardly ever needed. Hence, positive values of
%      SHIFTA/B indicate a shift to the rigth, and negative values a
%      shift to the left, whereas SHIFTDIM uses the opposite convention.
% 
%   2) SINGLETON EXPANSION (SX)
%      Whenever a dimension of either A or B is singleton and the
%      corresponding dimension of the other array is not, the mismatch is
%      fixed by virtually replicating (or diminishing to length 0) the
%      array along that dimension.
% 
%   BAXFUN applies elementwise operations. Matrix multiplications with AX
%   enabled can be performed using MULTIPROD (MATLAB Central, file #8773).
%   BAXFUN calls the builtin function BBXFUN, available in MATLAB R2007a.
%   For earlier MATLAB releases, use Schwarz's replacement of BSXFUN
%   (MATLAB Central, file #23005) 
%
%   Examples:
%       Subtracting the column means from a matrix A
%           a = magic(5); % ................. 5×5
%           a = baxfun(@minus, a, mean(a)); % 5×5
%
%       Subtracting the matrix means from 10 matrices contained in A
%           a = rand(3, 3, 10); % ............... 3×3×10
%           means = mean(reshape(a, 9, 10)); % .... 1×10
%           a = baxfun(@minus, a, means, 0, 1); % 3×3×10
% 
%       Multiplying matrix A by each element of B, i.e. multiplying each
%       element of A by each element of B (all possible combinations)
%           a = [1 2 3 4; 5 6 7 8]; % ....... 2×4
%           b = [1 10 100]; % ................. 1×3
%           c = baxfun(@times, a, b, 0, 1); % 2×4×3
%
%           c(:,:,1) =      1     2     3     4
%                           5     6     7     8
%
%           c(:,:,2) =     10    20    30    40
%                          50    60    70    80
%
%           c(:,:,3) =    100   200   300   400
%                         500   600   700   800
% 
%   See also  MULTIPROD (MATLAB Central, file #8773), BSXFUN, SHIFTDIM 
 
% $ Version: 1.2 $
% CODE      by:          Paolo de Leva
%                        (Univ. of Rome, Foro Italico, IT)      2009 Jan 31
% COMMENTS  by:          Code author                            2009 Feb 23
% OUTPUT    tested by:   Code author                            2009 Feb 23
% -------------------------------------------------------------------------

% Allow 3 to 5 input arguments
% This raises a warning R2016a and possibly errors in later versions.
%    error( nargchk(3, 5, nargin) );  


% Checking number of arguments
    switch nargin
        case 3
            c = bsxfun(f, a, b);
            return
        case 4
            shiftB = 0;
    end

% Checking for gross input errors
    if     1~=numel(shiftA) ||    1~=numel(shiftB) || ...
            ~isreal(shiftA) ||     ~isreal(shiftB) || ...
         ~isnumeric(shiftA) ||  ~isnumeric(shiftB) || ...
        shiftA~=fix(shiftA) || shiftB~=fix(shiftB) || ...
          ~isfinite(shiftA) ||   ~isfinite(shiftB)
        error('BAXFUN:InvalidShiftsize', ...
        'SHIFTA and SHIFTB must be finite integer scalars');
    end

% Processing
    if shiftA
        a = shiftdim(a, -shiftA);
    end
    if shiftB
        b = shiftdim(b, -shiftB);
    end
    c = bsxfun(f, a, b);
