% PCOLOR_CENTERED(first 3 basic pcolor arguments)
% 
% Makes pcolor actually plot the matrix you give it (e.g. give it a 10x10
% matrix, plot a 10x10 matrix, not a 9x9 matrix).

% function PC = pcolor_centered(X,Y,D)
function PC = pcolor_centered(varargin)

if nargin == 3

X = varargin{1};
Y = varargin{2};
D = varargin{3};

dx = X(1,2) - X(1,1);
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];

dy = Y(2,1) - Y(1,1);
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];

elseif nargin == 1

D = varargin{1};
[X,Y] = meshgrid(1:size(D,2), 1:size(D,1));

dx = 1;
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];

dy = 1;
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];

else
    error('1 or 3 inputs expected')
end

D_ = [D , nan(size(D,1),1) ; ...
      nan(1,size(D,2)) , nan ];

PC = pcolor(X_ - dx/2, Y_ - dy/2, D_);

shading flat

end
