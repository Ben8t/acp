#include <iostream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "function_acp.h"
#include <stdio.h> /* stderr */

#include "stack-c.h" /* Provide functions to access to the memory of Scilab */
#include "call_scilab.h" /* Provide functions to call Scilab engine */
#include "api_scilab.h"
using namespace std;

int main(){
    intro();
    string s;
    cout << "Entrer un fichier ('fichier.csv') : ";
    cin >> s;
    cout << endl;
    Matrice X(s);
    message("Matrice de base");
    X.display();
    message("Recentrage");
    X.refocus().display();
    message("Matrice de covariance");
    X.sigma().display();
    message("Plus grande valeur propre");
    cout << X.sigma().puissanceIter(0.05).val;
    message("Plus grand vecteur propre associe");
    displayVecteur(X.sigma().puissanceIter(0.05).vect);

    message("TEST DEFLATION");
    X.sigma().deflation().display();

    //REPRESENTATION DES INDIVIDUS

    Matrice ResPI = X.sigma().deflation();

    Matrice Y = ShowIndiv(X.refocus(), ResPI);
    message("Coord des individus");
    Y.display();

    Matrice C = ShowVar(X.refocus(), ResPI);
    message("Cercle de correlation");
    C.display();

    Matrice Contribution = Contrib(ResPI);
    message("Contribution");
    Contribution.display();

    bool c = choose("Voulez vous afficher les graphiques");
    if(c){
    ScilabExec(Y,C);
    }

    return 0;
}


//RESULT SIGMA

//12.469444	10.840278	2.973333	5.1125	0.1375
//10.840278	9.773611	4.592083	5.8625	0.075
//2.973333	4.592083	13.225	10.21875	0.4125
//5.1125	5.8625	10.21875	8.625	0.8125
//0.1375	0.075	0.4125	0.8125	9.75



//Characteristic Polynomial:
//  -x^4 + 40.6x^3 - 344.8x^2 + 3.753x + 1.1676
//
//Real Eigenvalues:   { -0.052854545539730624 ;  0.06415318495563642 ;  12.077546201161022 ;  28.511155159423076 }
//
//Eigenvectors:
//
// for Eigenvalue -0.052854545539730624:
//   [ -0.9140925178230374 ; 0.7941480166540801 ; -0.844852691468541 ; 1 ]
//
// for Eigenvalue 0.06415318495563642:
//   [ 0.6534306418440988 ; -1.0894521012887426 ; -0.5474117449991875 ; 1 ]
//
// for Eigenvalue 12.077546201161022:
//   [ -1.7388712441869858 ; -1.130952103964358 ; 2.0019416882426837 ; 1 ]
//
// for Eigenvalue 28.511155159423076:
//   [ 1.064005727303135 ; 1.04581632479926 ; 1.0154826924725278 ; 1 ]
//
//
//All tests OK


//Characteristic Polynomial:
//  -x^4 - 81.25999999999982x^3 + 1128.0031999999972x^2 - 44.45324099999959x + 0.31187198000088384
//
//Real Eigenvalues:   { -93.34884237523087 ;  0.00912911588808098 ;  0.03037214056863023 ;  12.049341118774343 }
//
//Eigenvectors:
//
// for Eigenvalue -93.34884237523087:
//   [ 1.0610314520927462 ; 1.046666926986179 ; 1.016297536009236 ; 1 ]
//
// for Eigenvalue 0.00912911588808098:
//   [ -2.245735472568279 ; 2.3819912953244446 ; -1.08997909259296 ; 1 ]
//
// for Eigenvalue 0.03037214056863023:
//   [ -0.10351342292943713 ; -0.18279893921650542 ; -0.6876777858025527 ; 1 ]
//
// for Eigenvalue 12.049341118774343:
//   [ -1.756279557238498 ; -1.1384540019380205 ; 2.022389053154941 ; 1 ]
//
//
//All tests OK
