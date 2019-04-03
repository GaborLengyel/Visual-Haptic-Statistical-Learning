function B = myfindfirst(A)
% B = FINDFIRST(A)
%
% Look for the row-indices of a first non-zero element(s) for all columns
% in the array.

B=zeros(1,size(A,1));
for j=1:size(A,2)
  Bj = find(A(:,j), 1, 'first');
  if ~isempty(Bj); B(j)=Bj; end
end


end
