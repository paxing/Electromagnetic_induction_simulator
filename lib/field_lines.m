% entry : B_field matrix
% return: vector lines plot

function hLines=field_lines(Bx_field,By_field)


FX=Bx_field;
FY=By_field;


[m,n]=size(Bx_field);

%logical vector (ignoring divergent point in gradient)
validColumns = all(isfinite(FX) & isfinite(FY));
[X,Y] = meshgrid(1:1:m,1:1:n);
%plot the vector lines
hLines = streamslice(X(:,validColumns),Y(:,validColumns),FX(:,validColumns),FY(:,validColumns));
set(hLines,'Color','k');
