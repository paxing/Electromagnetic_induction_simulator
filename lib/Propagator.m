%************************************************************
% Date : automne 2018
% auteur : Paul Xing

% entrée: 
% position_s : vecteur position source
% position_p : vecteur position point p
%
%
% sortie:
% P_tenseur: matrice 3x3, tenseur de propagation
%************************************************************



function P_tenseur=Propagator(position_s,position_p)

% calcul de la distance relative entre les 2 vecteurs
vecteur_R=(position_p-position_s);

%calcul de la norme du vecteur distance relative
norme_R=norm(vecteur_R);

%normalisation du vecteur de distance relative
vecteur_R=vecteur_R/norm(vecteur_R);
Id=eye(length(position_p));

%calcul du produit tensoriel
N_tenseur = 3*kron(vecteur_R,transpose(vecteur_R))-Id;

%calcul du propagateur (voir rapport)
P_tenseur=N_tenseur/(4*pi*norme_R^3);




