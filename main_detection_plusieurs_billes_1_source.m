clear all
clc
addpath("lib\")

%****************caractéristiques de la source****************************
%
N=100; %nombre de tour
I=1; % courant en ampère
A=pi*(0.01)^2; %surface de la bobine
ms=N*I*A; % moment dipolaire de la source
mu_0=4*pi*10^(-7);
omega=100000*2*pi; %fréquence angulaire
%**************************************************************************



%******************caractéristique de la bille*****************************
%
a=10/1000; %rayon en m
vs=4*pi*a^3/3; %volume de la bille
sigma=10^(6);
delta=sqrt(2/(mu_0*sigma*omega));

%calcul de la polarisabilité
alpha_0=-2*i/5*a^2*vs/delta^2;
%**************************************************************************



%******************initialisation du maillage******************************
%
M=301;
N=301;
x=linspace(0,0.3,M); %metre
y=linspace(0,0.3,N);%metre
%taille de la zone 30cm par 30cm max
%%

%*******************Positionnement des billes******************************

%position et orientation de la source
position_s=[0.150;0;0]; %en mètre
n_vec=[0;1;0]; %vecteur d'orientation



%calcul du moment dipolaire de la source
ms_vecteur=ms*n_vec;
bs_vecteur=mu_0*ms_vecteur;

%position des billes
N_bille=6;
position_1=[0.025;0.025;0]; % en mètre
position_2=[0.075;0.225;0];
position_3=[0.275;0.200;0];
position_4=[0.075;0.200;0];
position_5=[0.175;0.050;0];
position_6=[0.150;0.075;0];

%identation des positions dans un vecteur de matrices
position_obj=[position_1 position_2 position_3 position_4 position_5 position_6]; %en mètre
%polarisation isotrope identique pour toutes les billes
alpha=alpha_0*eye(3);
%identation des polarisation dans un vecteur de matrice
alpha=[alpha alpha alpha alpha alpha alpha];


%%
%***********************calcul du champ************************************
%calcul du champ en l'absence de bille (polarisation nulle)
[background_matrix,x_background,y_background]=B_field_multibille(x,y,position_s,bs_vecteur,position_obj,zeros(3,3*N_bille));

%calcul du champ en présence des billes
[space_matrix,x_matrix,y_matrix]=B_field_multibille(x,y,position_s,bs_vecteur,position_obj,alpha);

%%
%*****************************graphique************************************
%

wt=pi/4;

    x_field=real(x_matrix-x_background)*cos(wt)+imag(x_matrix-x_background)*sin(wt);
    y_field=real(y_matrix-y_background)*cos(wt)+imag(y_matrix-y_background)*sin(wt);
    field_norm=sqrt(x_field.^2+y_field.^2);
    field_image=imagesc(log10(field_norm));
hold on
c=colorbar;
colormap jet
c.Label.String = 'Log(B)';
xlabel('mm')
ylabel('mm')

hLinesplot=field_lines(imag(x_matrix),imag(y_matrix));
set(gca,'YDir','reverse');
xlim([0 N])
ylim([0 M])
set(gca,'FontSize',16)
x0=10;
y0=10;
width=900;
height=800;
set(gcf,'position',[x0,y0,width,height])
print -djpeg -r300 field_plusieurs_billes_final_diff
%%
figure(2)
hold on
field_image=imagesc(log10(space_matrix));
c=colorbar;
c.Label.String = 'Log(B)';
xlabel('mm')
ylabel('mm')
colormap jet
phase=pi/4;
hLinesplot=field_lines(real(x_matrix)*cos(phase)+imag(x_matrix)*sin(phase), real(y_matrix)*cos(phase)+imag(y_matrix)*sin(phase));
set(gca,'YDir','reverse');
xlim([0 M])
ylim([0 N])
set(gca,'FontSize',16)
x0=10;
y0=10;
width=900;
height=800;
set(gcf,'position',[x0,y0,width,height])
print -djpeg -r500 field_plusieurs_billes2
%%

%%
%*************************Animation en gif*****************************

% gif avec signal total
% Cette section permet d'enregister une animation de l'évolution
% temporelle du champ magnétique total
h = figure;

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for wt = 0:0.01*pi:2*pi
    
    x_field=real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt);
    y_field=real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt);
    field_norm=sqrt(x_field.^2+y_field.^2);
    field_image=imagesc(log10(field_norm));
    colormap jet
    c=colorbar;
    caxis([-9 1])
    c.Label.String = 'Log(B)';
    hold on
    field_lines(real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));
    
    xlabel('mm')
    ylabel('mm')
    xlim([0 M])
    ylim([0 N])
    set(gca,'FontSize',14)
    x0=10;
    y0=10;
    width=600;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'YDir','reverse');
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if wt == 0 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
      end 
end
  fprintf('GIF complet')


%%
%*************************Animation en gif*****************************

% gif avec signal total
% Cette section permet d'enregister une animation de l'évolution
% temporelle du champ magnétique total
h = figure;

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated_diff.gif';
for wt = 0:0.01*pi:2*pi
    
    x_field=real(x_matrix-x_background)*cos(wt)+imag(x_matrix-x_background)*sin(wt);
    y_field=real(y_matrix-y_background)*cos(wt)+imag(y_matrix-y_background)*sin(wt);
    field_norm=sqrt(x_field.^2+y_field.^2);
    field_image=imagesc(log10(field_norm));
    colormap jet
    c=colorbar;
    caxis([-9 1])
    c.Label.String = 'Log(B)';
    hold on
    field_lines(real(x_matrix-x_background)*cos(wt)+imag(x_matrix-x_background)*sin(wt), real(y_matrix-y_background)*cos(wt)+imag(y_matrix-y_background)*sin(wt));
    
    xlabel('mm')
    ylabel('mm')
    xlim([0 M])
    ylim([0 N])
    set(gca,'FontSize',14)
    x0=10;
    y0=10;
    width=600;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'YDir','reverse');
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if wt == 0 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
      end 
end
  fprintf('GIF complet')


