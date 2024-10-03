% omega = vorticity2D(U,V,dx,dy)
% 
% Calculates vorticity "omega" from a 2D current field {U,V} using some
% sort of central finite difference method.
% 
% IN:   U  = matrix of x-direction current
% IN:   V  = matrix of y-direction current (same size as U)
% IN:   dx = x step
% IN:   dy = y step
% 
% OUT:  omega = matrix (same size as U and V), vorticity

function omega = vorticity2D(U,V,dx,dy)

if (size(U,1) == size(V,1)) && (size(U,2) == size(V,2))
else
    error('U and V must be the same size.')
end

U_nan = nan(size(U,1) + 2, size(U,2) + 2);
U_nan(2:(end-1),2:(end-1)) = U;
U = U_nan;
clear U_nan
V_nan = nan(size(V,1) + 2, size(V,2) + 2);
V_nan(2:(end-1),2:(end-1)) = V;
V = V_nan;
clear V_nan

I = size(U,2);
J = size(U,1);
IJ = I*J;

% Make sparse matrices to implement this instead of slow loops.
SO_i = spdiags(repmat( [-1 0 1],I,1),[-1 0 1],speye(I)); % Sparse Operator matrix for i-dimension
SO_j = spdiags(repmat( [-1 0 1],J,1),[-1 0 1],speye(J)); % Sparse Operator matrix for j-dimension
A0_i = speye(I)*0;
A0_j = speye(J)*0;

% Build the big matrices (no preallocation because I don't want the headache of juggling this many indices):
SO_V = [];
for jj = 1:J
    SO_V = [SO_V; repmat(A0_i,1,jj-1), SO_i, repmat(A0_i,1,J-jj)];
end
SO_U = [];
for ii = 1:I
    SO_U = [SO_U; repmat(A0_j,1,ii-1), SO_j, repmat(A0_j,1,I-ii)];
end

U = reshape(U ,IJ,1);
V = reshape(V',IJ,1);

omega = reshape(SO_V*V,I,J)'/(2*dx) - reshape(SO_U*U,J,I)/(2*dy);
omega = omega(2:(end-1),2:(end-1));

end

% https://scicomp.stackexchange.com/questions/21915/discrete-definitions-of-curl-nabla-times-f
% https://en.wikipedia.org/wiki/Finite_difference_method

%% DEMONSTRATION

% close all
% 
% [X,Y] = meshgrid(-10:10,-10:10);
% 
% % U = sin(Y/5) + 5;
% % V = X/10;
% % U = sin(-Y/20).*exp(-0.05*X.^2 -0.05*Y.^2);
% % V = sin(X/20).*exp(-0.05*X.^2 -0.05*Y.^2);
% U = (Y+5).^2;
% V = (X-3).^2;
% 
% dx = 1; dy = 1;
% 
% OMEGA = vorticity2D(U,V,dx,dy);
% 
% figure
% 
% subplot(121)
% pcolor(X,Y,( sqrt( (U).^2 + (V).^2) )); hold on; shading flat
% quiver(X, Y, (U), (V),'color',[1 1 1]);
% 
% 
% subplot(122)
% pcolor(X,Y,OMEGA); shading flat
% colorbar
% 
% set(gcf,'Position',[255 378 1143 420])
% 
% 
% % % In Mathematica:
% % u = x y;(*or whatever*)
% % v = x y;
% % VectorPlot[{u, v}, {x, -10, 10}, {y, -10, 10}]
% % D[v, x] - D[u, y]
% % DensityPlot[%, {x, -10, 10}, {y, -10, 10}]
