%************************************************************
%date: mars 2019
%auteur: Paul Xing
%
%entrée:
%
% x:vecteur, maillage en x
% y: vecteur, maillage en y
%
% position_s: vecteur des vecteurs position de la source
%   (les vecteurs positions doivent être des vecteurs colonnes)
%
% bs_vecteur: vecteur des vecteurs du moment dipolaire de la source
%   (les vecteurs moment dipolaire doivent être des vecteurs colonnes)
%
% position_obj: vecteur des vecteurs positions des billes (matrice)
%   (les vecteurs positions doivent être des vecteurs colonnes)
%
% alpha_0 : vecteur des tenseurs de polarisabilité des billes
%
%
%
%sortie:
% space_matrix: matrice contenant la norme du champ magnétique
% x_matrix: matrice contenant la composante en x du champ
% y_matrix: matrice contenant la compoante en y du champ
%************************************************************





%************************************************************

function [space_matrix,x_matrix,y_matrix]= B_field_multibille_source(...
    x,y,position_s,bs_vecteur,position_obj,alpha_0)

%taille du maillage
M=length(x);
N=length(y);


%************************************************************
%initialisation des matrices
%************************************************************

%norme du champ magnétique selon la position
space_matrix=zeros(M,N);

%composante du champ selon x
x_matrix=zeros(M,N);

%composante du champ selon y
y_matrix=zeros(M,N);

%initialisation du vecteur de vecteurs qui contient la valeur du champ à la
%position de chaque bille causé par l'interaction des sources
B_s=[];

%détermine la dimension d du vecteur position et le nombre L de bille
[d,L]=size(position_obj);

%détermine la dimension h du vecteur position et le nombre K de source
[h,K]=size(position_s);

%création d'une matrice identité de taille d*L
vec1=speye(d*L,d*L);

%initialise la matrice de polarisation pour l'ensemble des billes
alpha=sparse(d*L,d*L);

%initialise la matrice P
P_N=sparse(d*L,d*L);
%************************************************************



%************************************************************
%Calcule de l'interaction des sources avec les billes
%************************************************************

for l=1:L %L est le nombre de billes
    
    %calcul du champ de chaque bille causé par l'interaction des sources
    B_sT=0;%initialisation du champ
    
    
    %somme la propagation de toutes les sources à chaque bille
    for s=1:K
        
        %calcul la contribution de la source s à k pour la bille l
        B_sn=Propagator(position_s(:,s),position_obj(:,l))*bs_vecteur(:,s);
       
        B_sT=B_sT+B_sn;%ajoute la valeur de la contribution de chaque bille
        %itérée
    end
    
    %compilation de la valeur du champ magnétique à la position de chaque
    %bille
    B_s=[B_s; B_sT];
    
    
    %matrice de polarisation des L billes
    alpha(1+d*(l-1):d+d*(l-1),1+d*(l-1):d+d*(l-1))= ...
        alpha_0(1:d,1+d*(l-1):d+d*(l-1));
    %************************************************************
    
    
    
    
    %************************************************************
    %calcul l'interaction entre chaque bille
    %************************************************************
    
    
    for k=1:L %L est le nombre total de bille
        if k~=l %permet d'exclure l'interaction d'une bille avec elle-même
           
            %calcul du propagateur entre les différentes billes
            P_N(1+d*(l-1):d+d*(l-1),1+d*(k-1):d+d*(k-1))= ...
                Propagator(position_obj(:,l),position_obj(:,k));   
        end
    end
end
%************************************************************




%************************************************************
%calcul de la matrice (voir rapport)
%************************************************************
A=vec1-transpose(alpha*P_N);

%calcul du champ à la position de chaque bille
%(résolution du système)
B_N=A\B_s;
%************************************************************



%************************************************************
%calcul du champ à chaque point du maillage
%************************************************************

for i=1:N % parcours les lignes (en y)
    for j=1:M %parcours les colonnes (en x)
        
        %propagateur des sources au pt de détection
        B_sd=0;
        %somme le champ de toutes les sources
        for s=1:K
            B_sdn=Propagator(position_s(:,s),[x(j);y(i);0])*bs_vecteur(:,s);
            B_sd=B_sd+B_sdn;
        end
        
        
        
        %propagateur de l'objet au pt de détection en tenant compte de la
        %polarisabilité
        
        aP_od=[]; %initialisation du vecteur de matrice
        
        for l=1:L
            %propagation de chaque bille au pt de détection
             aP_od= [aP_od alpha_0(1:d,1+d*(l-1):d+d*(l-1))*Propagator( ...
                 position_obj(:,l),[x(j);y(i);0])];
        end
        %valeur du champ causée par la présence des billes
        B_bille=aP_od*B_N;
        % comme on fait le produit scalaire entre un vecteur de matrice et
        % un vecteur de vecteurs, la somme se fait automatiquement
        
        
        %calcul du champ total
        if isfinite(B_bille)
            B_d=B_sd+B_bille; %permet de prendre en compte que les valeurs
            %définies
        else
             B_d=B_sd; %ignore les valeurs qui sont singulières
        end
        
        space_matrix(i,j)=norm(B_d); %calcul de la norme du champ
        %(valeur purement réelle)
        
        x_matrix(i,j)=B_d(1); %composante selon x du champ(valeur complexe)
        y_matrix(i,j)=B_d(2);% composante selon y dy champ(valeur complexe)
        
    end
end
