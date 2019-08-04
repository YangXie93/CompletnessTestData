#include <iostream>
#include <vector>
#include <list>
#include <Rcpp.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

//[[Rcpp::export]]

vector<int> randomContigs(int minContigLength,int meanContigLength,int covering){
    vector<int> res;
    default_random_engine generator;
    poisson_distribution<int> distribution(meanContigLength -minContigLength);
    int tmp;
    while(covering > 0){
        tmp = distribution(generator) +minContigLength;
        distribution.reset();
        if(tmp >= covering){
            res.push_back(covering);
            covering = 0;
        }
        else{
            covering -= tmp;
            res.push_back(tmp);
        }
    }
    return res;
}

//[[Rcpp::export]]

vector<int> randomSpaces(int numSp,int free){

    srand(time(NULL));
    vector<int > spaces;
    vector<int > nums;
    spaces.reserve(numSp);
    nums.reserve(numSp);
    int relativeSp = numSp/free;
    default_random_engine generator;
    poisson_distribution<int> distribution(relativeSp);
    for(int n = 0;n < numSp;n++){
        spaces.push_back(n);
        nums.push_back(n);
    }
    vector<int>::iterator it = spaces.begin();
    vector<int>::iterator n;
    while(numSp > 0){
        n = nums.begin();
        n += rand() % (int) nums.size();
        if((int) nums.size() > 1){
            (*(it +(*n))) = distribution(generator) % free;
            distribution.reset();
            free -= (*(it +(*n)));
            nums.erase(n);
        }
        else{
            (*(it +(*n))) = free;
        }
        numSp--;

    }
    return spaces;
}


//[[Rcpp::export]]

vector<int> fromWhichHowMany(int minContigLength,int totalLength,vector<int> lengths,double partCovered){
    default_random_engine generator;
    int needed = (int) (totalLength * partCovered);
    vector<int> res;
    res.reserve(lengths.size());
    int share;
    int tmp;
    for(vector<int>::iterator it = lengths.begin();it != lengths.end();it++){
        if((*it) >= minContigLength){
            share = round(needed *((*it)/(double)totalLength));
            poisson_distribution<int> distribution(share);
            tmp = distribution(generator) %needed;
            distribution.reset();
            needed -= tmp;
            res.push_back(tmp);
        }
        else{
            res.push_back(0);
        }

    }
    return res;
}


//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

Rcpp::List mkContigs(Rcpp::IntegerVector& lengths,int minContigLength,int meanContigLength, double partCovered = 0.8,bool sameContigLength = true){

    Rcpp::List res;
    int times;
    int rest;
    vector<int> covered;
    Rcpp::IntegerVector tmp1;
    Rcpp::IntegerVector tmp2;
    vector<int> spaces;
    vector<int> contigs;
    vector<int>::iterator co;
    vector<int>::iterator sp;
    bool swich;
    int bases;
    for(Rcpp::IntegerVector::iterator i = lengths.begin(); i != lengths.end();i++){
        bases = (int) ((*i) *partCovered);
        rest = (*i) - bases;                           // Problem: was wenn es nicht aufgeht?
        if(sameContigLength){
            times = bases/meanContigLength;
            vector<int> v(times,meanContigLength);
            contigs = v;
            rest += bases % meanContigLength;
        }
        else{
            contigs = randomContigs(minContigLength,meanContigLength,bases);
        }
        spaces = randomSpaces((int)contigs.size() +1,rest);
        swich = true;
        co = contigs.begin();
        sp = spaces.begin();
        tmp1 = Rcpp::IntegerVector::create();
        tmp2 = Rcpp::IntegerVector::create();
        for(int i = 0;i < (int) (contigs.size() << 1) +1;i++){
            if(swich){
                tmp1.push_back((*sp));
                tmp2.push_back(0);
                swich = false;
                sp++;
            }
            else{
                tmp1.push_back((*co));
                tmp2.push_back(1);
                swich = true;
                co++;
            }

        }
        res.push_back(tmp1);
        res.push_back(tmp2);
    }
    return res;
}
