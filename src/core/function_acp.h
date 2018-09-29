
//acp_function.h



#ifndef function_acp_h
#define function_acp_h


#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <stdio.h> /* stderr */

#include "stack-c.h" /* Provide functions to access to the memory of Scilab */
#include "call_scilab.h" /* Provide functions to call Scilab engine */
#include "api_scilab.h"

//Création d'un type qui permet de retourner un couple : valeur propre et vecteur propre associée.
struct ValVect{
    double val;
    std::vector<double> vect;
};


//Classe matrice : définition de nombreuses fonctions (manipulation de matrice, puissanceIterée, déflation, etc...
class Matrice{
    private :
        std::vector<std::vector<double> > m;
        int lon;
        int lar;

    public :

    //CONSTRUCTEUR
        //Constructeur à partir d'un fichier
        Matrice(std::string file){
            std::ifstream in(file.c_str());
            std::string a;
            int i(0);
            std::vector<std::string> col(0);
                if(in.fail()){
                    std::cerr << "Erreur avec le fichier : " << file;
                }
                else{
                    while(!in.eof()){
                        getline(in, a, '\n' );
                        col.push_back(a);
                        i++;

                    }
                    in.close();
                }
                for(int a(0);a<col.size();a++){
                    std::string s = col[a];
                    size_t pos = 0;
                    std::vector<double> vd;
                    double d = 0.0;

                    // convert ',' to ' '
                    while (pos < s.size ()){
                        if ((pos = s.find_first_of (';',pos)) != std::string::npos){
                            s[pos] = ' ';
                        }
                    }
                    std::stringstream ss(s);
                    while (ss >> d){//Convert & Push_back string to double
                        vd.push_back (d);
                    }

                    m.push_back(vd);
                    lar = vd.size();
                    lon = m.size();
                }
        };

        //Constructeur par défaut a 1
        Matrice(int c,int l){
            std::vector<std::vector<double> > A(c,std::vector<double>(l,1));
            m = A;
            lon = c;
            lar = l;
        }


    // DESTRUCTEUR

    //????//

    //METHODE

        //Afficher matrice ?? REMPLACER PAR COUT EN SURCHARGANT << ??
        void display() const;

        //Retourne le nombre de ligne
        int getLon();

        //Retourne le nombre de colonne
        int getLar();

        //Retourne le "tableau" (matrice)
        std::vector<std::vector<double> > getMatrice();

        //Modifier la matrice
        void modifMatrice(int i, int j, double c);

        //Calcul transposé
        Matrice transp();

        //Recentrage
        Matrice refocus();

        //Produit avec une matrice
        Matrice produit(Matrice B); // ?? SURCHARGE OPERATEUR * ??//

        //Sigma (Matrice de Variance-Covariance)
        Matrice sigma();

        //Produit d'une matrice avec un vecteur
        std::vector<double> produitVecteur(std::vector<double> x);

        //Puissance itérée,  On ne peut l'appeller que sur une matrice carrée ?? Faire une sous classe ou seulement gestion erreur??
        ValVect puissanceIter(double presicion);


        //Déflation : retourne un tableau contenant les vecteurs propres et valeurs propres associées
        // Application sur la matrice de Variance-Covariance
        Matrice deflation();

        Matrice produitVectTransp(std::vector<double> x, double lambda);

};


//Calcul la moyenne d'un vecteur
double mean(std::vector<double> x);

//Calcul de la variance
double variance(std::vector<double> x);

//Affichage vecteur
void displayVecteur(std::vector<double> x);

//Calcul la norme2 d'un vecteur
double norme2(std::vector<double> x);

//Soustraction entre deux vecteurs
std::vector<double> soustracVecteur(std::vector<double> x, std::vector<double> y);

//Affiche un message
void message(std::string s);

//Retourne une matrice composée des coordonnées des nouveaux individus
Matrice ShowIndiv(Matrice Xc, Matrice V);

//Retourne les coordonnée pour faire le cercle de corrélation
Matrice ShowVar(Matrice Xc, Matrice V);

//Tableau contenant les contribtuions de chaque axe à l'inertie totale (page 15)
Matrice Contrib(Matrice V);

//Création d'une matrice à la main
Matrice CreaMain();

//Affichage des graphiques (individus et cercle de correlation ) avec l'API Scilab
void ScilabExec(Matrice I, Matrice V);

void intro();

bool choose(std::string s);
#endif // function_acp_h



