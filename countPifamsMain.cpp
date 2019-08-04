#include "countPifams.cpp"

using namespace std;

int main(){
    list<vector<int> > a;
    list<vector<int> > b;

    vector<int> va = {10,15,15,10,10};
    vector<int> vb = {-1,3,-1,5,-1};

    vector<int> xa = {20,50};
    vector<int> xb = {50,60};
    vector<int> c = {0};

    b.push_back(c);
    b.push_back(xa);
    b.push_back(xb);


    list<list<vector<int> > > pi;
    list<list<vector<int> > > con;

    for(int i = 0;i < 2;i++){
        a.push_back(va);
        a.push_back(vb);
    }
    con.push_back(b);
    con.push_back(b);
    pi.push_back(a);

    Rcpp::List res = countPifams(pi,con);

    return 0;
}
