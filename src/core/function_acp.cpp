//function_acp.cpp

#include "function_acp.h"
using namespace std;


//MATRICE METHODES
void Matrice::display() const{
    cout << endl;
    cout << "//DATA DISPLAY//" <<endl<<endl;
    for(int i(0);i<lon;i++){
        for(int j(0);j<lar;j++){
            cout << m[i][j]  << "\t";
        }
        cout << endl;
    }
}

int Matrice::getLon(){
    return lon;
}

int Matrice::getLar(){
    return lar;
}

vector<vector<double> > Matrice::getMatrice(){
    return m;
}

void Matrice::modifMatrice(int i, int j, double c){
    m[i][j]=c;
}

Matrice Matrice::transp(){
    Matrice T(lar,lon);
    for(int i(0);i<lar;i++){
        for(int j(0);j<lon;j++){
            T.m[i][j]=m[j][i];
        }
    }
    return T;
}

Matrice Matrice::refocus(){
    vector<double> G;
    Matrice R(lon,lar);
    for(int i(0);i<lar;i++){
        G.push_back(mean(this->transp().m[i]));
    }
    for(int i(0);i<lon;i++){
        for(int j(0);j<lar;j++){
            R.m[i][j]=m[i][j]-G[j];
        }
    }
    return R;
}

Matrice Matrice::produit(Matrice B){
    double temp(0);
    Matrice P(lon,B.lar);
    for(int i(0);i<lon;i++){
        for(int j(0);j<lon;j++){
            for(int k(0);k<B.lon;k++){
                temp = m[i][k]*B.m[k][j];
                P.m[i][j]=P.m[i][j]+temp;
            }
        }
    }
    return P;
}

Matrice Matrice::sigma(){
    Matrice Xc = this->refocus();
    Matrice XcT = this->refocus().transp();
    Matrice P = XcT.produit(Xc);
    Matrice S(P.lon,P.lon);
    for(int i(0);i<P.lon;i++){
        for(int j(0);j<P.lon;j++){
                S.m[i][j]=P.m[i][j]/lon;
        }
    }
    return S;
}

vector<double> Matrice::produitVecteur(vector<double> x){//Gestion erreur ??
    vector<double> res(x.size(),0);
    for(int i(0);i<lon;i++){
        for(int j(0);j<lon;j++){
            res[i] += m[i][j]*x[j];
        }
    }
    return res;
}

ValVect Matrice::puissanceIter(double eps){
    ValVect res;
    double erreur = eps+1;
    vector<double> xk(lon,1/sqrt(double(lon)));
    vector<double> uk(lon,1/sqrt(double(lon)));
    vector<double> yk(lon,1);
    double lambdak;
    int cmp=0;
    while(erreur>eps && cmp<1000){
        lambdak=0;
        //Xk+1 = A*xk
        uk=this->produitVecteur(xk);
        //yk = xk/||xk||
        for(int i(0);i<lon;i++){
                yk[i]=uk[i]/norme2(uk);
        }
        //lambda = yk'*A*yk
        for (int i(0);i<lon;i++){
            lambdak+=yk[i]*this->produitVecteur(yk)[i];
        }
        xk=yk;
        cmp++;
        //IL FAUT MODIFIE LE CALCUL DE L'ERREUR !!!!!!!!
        erreur -= 0.05;
    }
    res.val = lambdak;
    res.vect = xk;
    return res;
}

Matrice Matrice::produitVectTransp(vector<double> x, double lambda){
    Matrice R(x.size(),x.size());
    for(int i(0);i<x.size();i++){
        for(int j(0);j<x.size();j++){
            R.m[i][j]=x[i]*x[j]*lambda;
        }
    }
    return R;
}

Matrice Matrice::deflation(){
    Matrice R(lon,lon+1);
    Matrice B(lon,lon);
    Matrice TMP(lon,lon);
    TMP.m = this->m;
    double lambda = this->puissanceIter(0.05).val;
    vector<double> vp = this->puissanceIter(0.05).vect;
    //Déflation
    for(int w(0);w<lon;w++){
        for(int i=0;i<lon;i++){
            for(int j=0;j<lon;j++){
                //Stockage
                R.m[w][j]=vp[j];
                //B=A - lambda1 * transp(u1)*u1
                B.m[i][j] = TMP.m[i][j] - produitVectTransp(vp,lambda).m[i][j];
            }
            //Stockage
            R.m[w][lon]=lambda;
        }
    lambda = B.puissanceIter(0.05).val;
    vp = B.puissanceIter(0.05).vect;
    TMP = B;
    }
    //Les resultat son stocké dan sun tableau, les  colonnes 0:n-1 corespondent aux composantes des vecteurs propre (normalisé) et la colonne n contient les valeurs propre associé
    return R;
}


//END METHODE MATRICE


double mean(std::vector<double> x){
    double sum(0);
    for(int i(0);i<x.size();i++){
        sum+=x[i];
    }
    return sum/x.size();
}

double variance(std::vector<double> x){
    double sum(0);
    for(int i(0);i<x.size();i++){
        sum+=(x[i]-mean(x))*(x[i]-mean(x));
    }
    return sum/x.size();

}

void displayVecteur(vector<double> x){
    cout << endl;
    for(int i(0);i<x.size();i++){
        cout << x[i] << ' ';
    }
    cout << endl;
}

double norme2(std::vector<double> x){
    double sum(0);
    for(int i(0);i<x.size();i++){
            sum += x[i]*x[i];
    }
    return sqrt(sum);
}

vector<double> soustracVecteur(vector<double> x, vector<double> y){
    vector<double> s(x.size(),0);
    for(int i(0);i<x.size();i++){
        s[i]=x[i]-y[i];
    }
    return s;
}

void message(std::string s){
    cout << "\n\n\n\n\n\n\t\t////////////////////////////////// \t" << s << "\t //////////////////////////////////\n\n";
}

Matrice ShowIndiv(Matrice Xc, Matrice V){
    //Xc matrice rencentré et V matrice qui contient les valeurs et vecteurs propre.
    Matrice Y(Xc.getLon(), V.getLon());
    double c=0;
        for(int i(0);i<Xc.getLon();i++){
            for (int j(0); j<V.getLon(); j++){
                for(int k(0);k<V.getLon();k++){
                    c+=V.getMatrice()[j][k]*Xc.getMatrice()[i][k];
                }
            Y.modifMatrice(i,j,c);
            c=0;
            }
        }
    return Y;
}


Matrice ShowVar(Matrice Xc, Matrice V){
    Matrice C(Xc.getLar(), V.getLon());
    double tmp;
    for(int i(0);i<Xc.getLar();i++){
        for(int j(0);j<V.getLon();j++){
            tmp=sqrt(V.getMatrice()[i][V.getLon()])*V.getMatrice()[i][j]/sqrt(variance(Xc.transp().getMatrice()[i]));
            //cout << endl << V.getMatrice()[i][j];
            C.modifMatrice(j,i,tmp);
        }
        cout << endl;
    }
    return C;
}


Matrice Contrib(Matrice V){
    Matrice R(V.getLon(),3);
    double sum(0), sum2(0);

    for (int k(0);k<V.getLon();k++){
        sum+=V.getMatrice()[k][V.getLon()];
    }
    for(int i(0);i<V.getLon();i++){
        for (int j(0);j<3;j++){
            switch (j){
            case 0:
                R.modifMatrice(i,j,V.getMatrice()[i][V.getLon()]);
                break;
            case 1:
                R.modifMatrice(i,j,V.getMatrice()[i][V.getLon()]/sum);
                sum2+=V.getMatrice()[i][V.getLon()]/sum;
                break;
            case 2:
                R.modifMatrice(i,j,sum2);
                break;
            }
        }
    }
    return R;
}

Matrice CreaMain(){
    int lon, lar;
    double a;
    cout << endl << endl << "INPUT MATRICE" << endl << endl;
    cout <<"Entrer le nombre de ligne : "; cin >> lon;
    cout << "Entrer le noombre de colonne : "; cin >> lar;
    Matrice A(lon,lar);
    for(int i(0);i<A.getLon();i++){
        for(int j(0);j<A.getLar();j++){
            cout << "Valeur pour " << i << "-" << j << " : "; cin >> a;
            A.modifMatrice(i,j,a);
        }
    }
    A.display();
    return A;
}

void ScilabExec(Matrice I, Matrice V){

    double * i1 =new double[I.getLon()];
    double * i2 =new double[I.getLon()];
    double * v1 =new double[V.getLon()];
    double * v2 =new double[V.getLon()];
    for(int i(0);i<I.getLon();i++){
        i1[i]=I.getMatrice()[i][0];
        i2[i]=I.getMatrice()[i][1];
    }

    for(int i(0);i<V.getLon();i++){
        v1[i]=V.getMatrice()[i][0];
        v2[i]=V.getMatrice()[i][1];
    }
    // A simple call_scilab example


/****** INITIALIZATION **********/
#ifdef _MSC_VER
 if ( StartScilab(NULL,NULL,0) == FALSE )
#else
 if ( StartScilab(getenv("SCI"),NULL,NULL) == FALSE )
#endif
  {
   fprintf(stderr,"Error while calling StartScilab\n");
  }

/****** ACTUAL Scilab TASKS *******/


char variableNamei1[] = "i1";
char variableNamei2[] = "i2";
char variableNamev1[] = "v1";
char variableNamev2[] = "v2";

createNamedMatrixOfDouble(pvApiCtx,variableNamei1,1,I.getLon(),i1);
createNamedMatrixOfDouble(pvApiCtx,variableNamei2,1,I.getLon(),i2);

createNamedMatrixOfDouble(pvApiCtx,variableNamev1,1,V.getLon(),v1);
createNamedMatrixOfDouble(pvApiCtx,variableNamev2,1,V.getLon(),v2);

SendScilabJob("plot2d(i1,i2,style=-2, axesflag=5);xtitle('Représentation des individus');");
system("pause");
SendScilabJob("plot2d(v1,v2,style=-2, frameflag=1, rect=[-1,-1,1,1], axesflag=5);xtitle('Cercle de corrélation');R=1;a=1;square(-a, -a, a, a);t = linspace(0, 2*%pi, 1001);x=R*cos(t);y=R*sin(t);plot2d(x,y)");
system("pause");
/****** TERMINATION **********/
 if ( TerminateScilab(NULL) == FALSE ) {
  fprintf(stderr,"Error while calling TerminateScilab\n");
 }
}

void intro(){
    cout << endl << endl;
    cout <<"\t\t\t\t\t\t" << "PROJET DATA SCIENCE - EISTI 2015/2016" << "\n\n";
    cout << "\t\t\t" << "AHAMMAD Ayoub, AUDRAIN Camille, HUCHE Mathieu, MULLEY Julien, PIMPAUD Benoit, SIMONI Vincent" << "\n\n\n";

    cout <<"\t\t\t\t\t\t" << "Analyse en Composantes Principales" << "\n\n\n\n";
}


bool choose(string s){
    char c;
    cout << "\n\n"<< s << " (y/n): ";
    cin >> c;
    if(c=='y'){
        return true;
    }
    else if(c=='n'){
        return false;
    }
    else{
        choose(s);
    }

}
