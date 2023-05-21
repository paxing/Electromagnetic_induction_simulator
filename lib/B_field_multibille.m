%date: mars 2019
%auteur: Paul Xing
%entr�e
% x:vecteur, maillage en x
% y: vecteur, maillage en y
% position_s: vecteur position de la source
% bs_vecteur: vecteur du moment dipolaire de la source
% position_obj: vecteur des vecteurs positions des billes (matrice)
% alpha_0 : vecteur des tenseurs de polarisabilit� des billes


function [space_matrix,x_matrix,y_matrix]=B_field_multibille(x,y,position_s,bs_vecteur,position_obj,alpha_0)

%taille du maillage
M=length(x);
N=length(y);


%initialisation des matrices
space_matrix=zeros(M,N);
x_matrix=zeros(M,N);
y_matrix=zeros(M,N);
B_s=[];
% d�termine la dimension du vecteur et le nombre d'objet L
[d,L]=size(position_obj);
vec1=speye(d*L,d*L);

%initialise la matrice de polarisabilit� contenant toutes les billes
alpha=sparse(d*L,d*L);

P_N=sparse(d*L,d*L);


for l=1:L %L est le nombre de bille
    
    %calcul du champ de chaque bille caus� par l'interaction de la source
    
    %selectionne chaque bille l un par un et calcule le propagateur et le
    %champ
    %compile das un vecteur de vecteur B_s le r�sultat pour chaque bille
    B_s=[B_s; Propagator(position_s,position_obj(:,l))*bs_vecteur];
    
    %vec1(1+d*(l-1):d+d*(l-1),1+d*(l-1):d+d*(l-1))=eye(d);
    
    %matrice de polarisation des L billes
    alpha(1+d*(l-1):d+d*(l-1),1+d*(l-1):d+d*(l-1))=alpha_0(1:d,1+d*(l-1):d+d*(l-1));
    
    %calcul l'interaction entre chaque bille
    for k=1:L
        if k~=l
            %calcul du propagateur entre les diff�rentes billes
            P_N(1+d*(l-1):d+d*(l-1),1+d*(k-1):d+d*(k-1))=Propagator(position_obj(:,l),position_obj(:,k));   
        end
    end
end

A=vec1-transpose(alpha*P_N);

%calcul du champ � la position de chaque bille
B_N=A\B_s;


%calcul du champ � chaque point du maillage
for i=1:M % parcours les lignes (en y)
    for j=1:N %parcours les colonnes (en x)
        
        %propagateur de la source au pt de d�tection
        P_sd=Propagator(position_s,[x(j);y(i);0]);
        
        %propagateur de l'objet au pt de d�tection
        aP_od=[];
        for l=1:L
            %propagation de chaque bille au pt de d�tection
             aP_od= [aP_od alpha_0(1:d,1+d*(l-1):d+d*(l-1))*Propagator(position_obj(:,l),[x(j);y(i);0])];
        end
        %valeur du champ caus�e par la pr�sence des billes
        B_bille=aP_od*B_N;
        
        %calcul du champ total
        
        if isfinite(B_bille) %permet de prendre en compte que les valeurs
            %d�finies
            B_d=P_sd*bs_vecteur+B_bille;
        else
             B_d=P_sd*bs_vecteur; %ignore les valeurs qui sont singuli�res
        end
        
        space_matrix(i,j)=norm(B_d); %calcul de la norme du champ
        x_matrix(i,j)=B_d(1); %composante selon x du champ
        y_matrix(i,j)=B_d(2);% composante selon y dy champ
        
    end
end
