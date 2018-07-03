function diagonal = diagProduct(A,B)
% diagonal = diagProduct(A,B)
% Efficient computation of diag(A*B) when A*B is square, 
% ie A is DxM, B is MxD.

diagonal = sum(A.*B',2);