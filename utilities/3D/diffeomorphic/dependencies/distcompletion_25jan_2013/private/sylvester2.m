function omega=sylvester2(sym_mat,asym_mat)
% omega=sylvester2(sym_mat,asym_mat)
% This function solves the Sylvester equation: Omega*sym_mat+sym_mat*Omega=asym_mat
% 
% Reference: 
% M. Journée, F. Bach, P.-A. Absil and R. Sepulchre, Low-rank optimization for semidefinite convex problems, arXiv:0807.4423v1, 2008
%

p=size(sym_mat,2);
[V,D]=eig(sym_mat);
O=zeros(p,p);
O2=V'*asym_mat*V;
for i=1:p-1,
    for j=i+1:p,
        A=D(i,i)+D(j,j);
        if A~=0,
            O(i,j)=O2(i,j)/A;
        else
           O(i,j)=0;
        end
    end
end

O=O-O';
omega=V*O*V';
