function [angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,x_matrix,y_matrix)


x_center=floor(M/2);
y_center=floor(N/2);
r = floor(M/2);


n_pts=900;
angle_vector=linspace(0, 2*pi, n_pts);
%crée un vecteur avec la position des pixels sur le cercle de rayon r
[x_polar,y_polar] = pol2cart(angle_vector, r);
x_polar=ceil(x_polar+x_center);
y_polar=ceil(y_polar+y_center);
B_field_radius=zeros(1,n_pts);

for i=1:n_pts
    angle=angle_vector(i);
    x_coor=x_polar(n_pts-i+1);
    y_coor=y_polar(n_pts-i+1);

     B_x=(x_matrix(y_coor,x_coor))*cos(angle);
     B_y=(y_matrix(y_coor,x_coor))*sin(angle);
    B_field_radius(i)=B_x+B_y;
end

