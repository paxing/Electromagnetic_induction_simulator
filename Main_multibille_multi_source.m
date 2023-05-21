clear all
clc
%%
addpath("lib\")
%****************caractéristiques de la source****************************
%
% On considère ici que toutes les sources ont les mêmes caractéristiques
% il est possible de considérer des sources différentes en définissant des
% variables N_1, N_2..., I_1,I_2....etc


N=100; %nombre de tour
I=1; % courant en ampère
A=pi*(0.01)^2; %surface de la bobine
ms=N*I*A; % moment dipolaire de la source
mu_0=4*pi*10^(-7);
omega=100000*2*pi; %fréquence angulaire
%**************************************************************************



%******************caractéristique de la bille*****************************
%On considère ici que toutes les billes ont les mêmes caractéristiques

a=10/1000; %rayon en m
vs=4*pi*a^3/3; %volume de la bille
sigma=10^(6);%conductivité
delta=sqrt(2/(mu_0*sigma*omega));%distance de pénétration

%calcul de la polarisabilité
alpha_0=-2*i/5*a^2*vs/delta^2;
%on définit ici une polarisabilité purement complexe qui correspond à une
%différence de phase


%**************************************************************************



%******************initialisation du maillage******************************
% Un maillage carré a été choisi ici

M=301;
N=301;
x=linspace(0,0.3,M); %metre
y=linspace(0,0.3,N);%metre
%la taille de la zone d'intérêt a été fixé ici à 30cm par 30cm. Il s'agit
%d'une taille réaliste qui peut être retrouvé en laboratoire


%%

%*******************Positionnement des sources*****************************

%nombre de sources, à mettre à jour si l'utilisateur choisi d'enlever ou
%d'ajouter des sources
N_source=4;

%position et orientation des sources
position_s1=[0.150;0;0]; %en mètre
position_s2=[0.150;0.300;0]; %en mètre
position_s3=[0.0;0.150;0]; %en mètre
position_s4=[0.300;0.150;0]; %en mètre
%note1: Attention à la définition des axes dans MATLAB

% l,utilisateur peut ajouter des sources en définissant position_s5..etc.
% les positions doivent respecter la taille de maillage

%vecteur d'orientation de chaque source avec PHASE
n_vec1=[0;i;0]; 
n_vec2=[0;-i;0]; 
n_vec3=[i;0;0]; 
n_vec4=[-i;0;0]; 
%note1: les vecteurs doivent être NORMALISÉE.
%note2: Attention à la définition des axes dans MATLAB
%note2: une orientation purement imaginaire i signifie un déphasage de pi/2
% par rapport à une orientation purement rélle de 1. Une différence de
% signe signifie une différence de phase de pi

 

%calcul du moment dipolaire des sources (sources identiques)
bs_vecteur1=mu_0*ms*n_vec1;
bs_vecteur2=mu_0*ms*n_vec2;
bs_vecteur3=mu_0*ms*n_vec3;
bs_vecteur4=mu_0*ms*n_vec4;

%identation des positions dans un vecteur de vecteurs
position_s=[position_s1 position_s2 position_s3 position_s4];
%NOTE: Ne pas oublier de rajouter les sources supplémentaires le cas
%échéant


%identation des moments dipolaires dans un vecteur de vecteur
bs_vecteur=[bs_vecteur1 bs_vecteur2 bs_vecteur3 bs_vecteur4];
%NOTE: Ne pas oublier de rajouter les sources supplémentaires le cas
%échéant

%%
%*******************Positionnement des billes******************************
%nombre de billes à mettre à jour si l'utilisateur choisi d'enlever ou
%d'ajouter des billes
N_bille=8;

%position des billes
position_1=[0.100;0.100;0]; % en mètre
position_2=[0.100;0.200;0];
position_3=[0.200;0.100;0];
position_4=[0.200;0.200;0];
position_5=[0.150;0.075;0];
position_6=[0.150;0.225;0];
position_7=[0.075;0.150;0];
position_8=[0.225;0.150;0];

%identation des positions dans un vecteur de vecteurs
position_obj=[position_1 position_2 position_3 position_4 position_5 position_6 position_7 position_8]; %en mètre
%NOTE: Ne pas oublier de rajouter les billes supplémentaires le cas
%échéant


%polarisation isotrope identique pour toutes les billes
%alpha=alpha_0*eye(3);
%NOTE: la notation tensorielle est utilisée ici

alpha_0=alpha_0*(1+i)/sqrt(2);
%polarisation pour une tige ellipsoide
alpha=[1*alpha_0  3*alpha_0 0;  3*alpha_0 2*alpha_0 0;0 0 0];

%identation des polarisation dans un vecteur de matrice (billes identiques)
alpha=[alpha alpha alpha alpha alpha alpha alpha alpha];
%NOTE: Ne pas oublier de rajouter les sources supplémentaires le cas
%échéan

%%
%***********************calcul du champ************************************
% le calcul ici considère l'ensemble des sources


%calcul du champ en l'absence de bille (polarisation nulle)
[background_matrix,x_background,y_background]=B_field_multibille_source(x,y,position_s,bs_vecteur,position_obj,zeros(3,3*N_bille));

%calcul du champ en présence des billes
[space_matrix,x_matrix,y_matrix]=B_field_multibille_source(x,y,position_s,bs_vecteur,position_obj,alpha);

%%
%***********************calcul du champ************************************
% calcul toutes les combinaisons possibles des sources

% les cellules ont été utilisé afin d'enregistrer les différentes
% combinaisons possible car elle permettent facilement d'enregistrer et
% d'appeler des matrices complètes

%ATTENTION: calcul très long...

%initialisation des cellules pour le champ total (norme et composantes x,y)
B_combnk=cell(linspace(N_source,N_source,N_source));
B_combnk_x=cell(linspace(N_source,N_source,N_source));
B_combnk_y=cell(linspace(N_source,N_source,N_source));


%initialisation des cellules pour le champ en l'absence de bille
B_combnk_bckgrd=cell(linspace(N_source,N_source,N_source));
B_combnk_x_bckgrd=cell(linspace(N_source,N_source,N_source));
B_combnk_y_bckgrd=cell(linspace(N_source,N_source,N_source));




%calcul itératif de toutes les combinaisons de sources possibles
for m=1:N_source
    %donne les combinaisons de sources possibles
    comb=nchoosek(1:N_source,m);
    [k,l]=size(comb); %k nombre de possibilités pour m sources
    indx=ones(1,N_source);
    
    for s=1:k
        h=comb(s,:);
        indx(1:l)=h;
        indxcell = num2cell(indx);
        
       
       [space_matrix,x_matrix,y_matrix]=B_field_multibille_source(x,y,position_s(:,h),bs_vecteur(:,h),position_obj,alpha);
       B_combnk{indxcell{:}}=space_matrix;
       B_combnk_x{indxcell{:}}=x_matrix;
       B_combnk_y{indxcell{:}}=y_matrix;
       
       [background_matrix,x_background,y_background]=B_field_multibille_source(x,y,position_s(:,h),bs_vecteur(:,h),position_obj,zeros(3,3*N_bille));
       B_combnk_bckgrd{indxcell{:}}=background_matrix;
       B_combnk_x_bckgrd{indxcell{:}}=x_background;
       B_combnk_y_bckgrd{indxcell{:}}=y_background;
    end
    
end


%%
%*****************données des combinaisons**********************
%
% cette section permet de séléctionner une combinaison spécifique de source
% les soruces doit être sélection dans l'ordre croissant.
% l'exemple ci-dessus permet de regarder la combinaison entre la source 1
% et 3.

space_matrix=B_combnk{1,2,3,4};
x_matrix=B_combnk_x{1,2,3,4};
y_matrix=B_combnk_y{1,2,3,4};

background_matrix=B_combnk_bckgrd{1,2,3,4};
x_background=B_combnk_x_bckgrd{1,2,3,4};
y_background=B_combnk_y_bckgrd{1,2,3,4};

%%
%*************************zone circulaire diff*****************************

%pour tracer une image, il faut sélection un argument de ««phase»»  wt car
%l'amplitude du signal varie en fonction du temps.
figure()
wt=pi/4; 



x_field=real(x_matrix-x_background)*cos(wt)+imag(x_matrix-x_background)*sin(wt);
y_field=real(y_matrix-y_background)*cos(wt)+imag(y_matrix-y_background)*sin(wt);
field_norm=sqrt(x_field.^2+y_field.^2);
    
[angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,abs(x_field),abs(y_field));
subplot(2,1,1);

field_image=imagesc(log10(field_norm));

hold on
c=colorbar;
caxis([-8.5 1])
colormap jet
c.Label.String = 'Log(B)';
xlabel('mm')
ylabel('mm')

hLinesplot=field_lines(x_field, y_field);
set(gca,'YDir','reverse');
xlim([0 M])
ylim([0 N])
set(gca,'FontSize',14)
x0=10;
y0=10;
width=400;
height=750;
set(gcf,'position',[x0,y0,width,height])

hold on
plot(x_polar,y_polar,'k','LineWidth',2)

subplot(2,1,2);
semilogy(angle_vector,abs(B_field_radius),'b','LineWidth',2)
ylabel('B (T)')
 set(gca,'XTick',0:pi/2:2*pi) 
 set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
 set(gca,'FontSize',14)
xlim([0 2*pi])
print -djpeg -r300 field_objet_aleatoire_polarisation_mixte

%%
%*************************zone circulaire total****************************

figure()
wt=pi/4;

[angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));
subplot(2,1,1);
field_image=imagesc(log10(space_matrix));
hold on
c=colorbar;
caxis([-9 1]);
colormap jet
c.Label.String = 'Log(B)';
xlabel('mm')
ylabel('mm')

hLinesplot=field_lines(real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));

set(gca,'YDir','reverse');
xlim([0 M])
ylim([0 N])
set(gca,'FontSize',14)
x0=10;
y0=10;
width=400;
height=750;
set(gcf,'position',[x0,y0,width,height])


hold on
plot(x_polar,y_polar,'LineWidth',2)

subplot(2,1,2);
plot(angle_vector,(B_field_radius)*10^9,'b','LineWidth',2)
ylabel('B (nT)')
 set(gca,'XTick',0:pi/2:2*pi) 
 set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
 set(gca,'FontSize',14)
xlim([0 2*pi])
print -djpeg -r300 field_circulaire_total_filiforme

%%
wt=pi/4;
figure()
[angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));
field_image=imagesc(log10(space_matrix));
hold on
c=colorbar;
caxis([-9 1]);
colormap jet
c.Label.String = 'Log(B)';
xlabel('mm')
ylabel('mm')

hLinesplot=field_lines(real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));

set(gca,'YDir','reverse');
xlim([0 M])
ylim([0 N])
set(gca,'FontSize',14)
x0=10;
y0=10;
width=900;
height=800;
set(gcf,'position',[x0,y0,width,height])
print -djpeg -r300 field_total1-2-3-4
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
  
  % animation avec évolution du signal autour de la zone circulaire
h = figure;

axis tight manual 
filename = 'testAnimated-test-phase_rot_2.gif';
n=0;
for wt = 0:0.01*pi:2*pi
    
    subplot(1,2,1);
    x_field=real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt);
    y_field=real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt);
    field_norm=sqrt(x_field.^2+y_field.^2);
    [angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));

    field_image=imagesc(log10(field_norm));
    colormap jet
    c=colorbar;
    caxis([-9 1])
    c.Label.String = 'Log(B)';
    hold on
    field_lines(real(x_matrix)*cos(wt)+imag(x_matrix)*sin(wt), real(y_matrix)*cos(wt)+imag(y_matrix)*sin(wt));
    
    hold on
    plot(x_polar,y_polar,'LineWidth',2,'Color','k')
    
    xlabel('mm')
    ylabel('mm')
    xlim([0 M])
    ylim([0 N])
    set(gca,'FontSize',14)
    x0=10;
    y0=10;
    width=1350;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'YDir','reverse');
    
    
    subplot(1,2,2);
    plot(angle_vector,log10(abs(B_field_radius)),'b','LineWidth',2)
    ylabel('Log(B)')
    ylim([-9 1])
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'FontSize',14)
    
    legend(strcat('\omega t= ', num2str(0.01*n,'%0.2f'),'\pi'))
    xlim([0 2*pi])
    n=n+1;
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

 % animation avec soustraction du signal des sources et avec signal autour
 % de la zone circulaire 
h = figure;

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated-test-forme-aleatoirev2.gif';
n=0;
for wt = 0:0.01*pi:2*pi
    
    subplot(1,2,1);
    x_field=real(x_matrix-x_background)*cos(wt)+imag(x_matrix-x_background)*sin(wt);
    y_field=real(y_matrix-y_background)*cos(wt)+imag(y_matrix-y_background)*sin(wt);
    field_norm=sqrt(x_field.^2+y_field.^2);
    [angle_vector,x_polar,y_polar,B_field_radius]=B_field_r(M,N,x_field,y_field);

    field_image=imagesc(log10(field_norm));
    colormap jet
    c=colorbar;
    caxis([-9 1])
    c.Label.String = 'Log(B)';
    hold on
    field_lines(x_field,y_field);
    
    hold on
    plot(x_polar,y_polar,'LineWidth',2,'Color','k')
    
    xlabel('mm')
    ylabel('mm')
    xlim([0 M])
    ylim([0 N])
    set(gca,'FontSize',14)
    x0=10;
    y0=10;
    width=1350;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'YDir','reverse');
    
    
    subplot(1,2,2);
    plot(angle_vector,log10(abs(B_field_radius)),'b','LineWidth',2)
    ylabel('Log(B)')
    ylim([-9 1])
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'FontSize',14)
    
    legend(strcat('\omega t= ', num2str(0.01*n,'%0.2f'),'\pi'))
    xlim([0 2*pi])
    n=n+1;
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