#include <iostream>
#include <vector>
#include <list>
//#include <Rcpp.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <algorithm>

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
    vector<int > tmp;
    spaces.reserve(numSp);
    tmp.reserve(numSp);
    int relativeSp = free/numSp;
    default_random_engine generator;
    poisson_distribution<int> distribution(relativeSp);
    int n;
    while(numSp > 0){
        if(numSp > 1){
            n = distribution(generator) % free;
            tmp.push_back(n);
            distribution.reset();
            free -= n;
        }
        else{
            tmp.push_back(free);
        }
        numSp--;

    }
    vector<int>::iterator it;
    while((int) tmp.size() > 0){
         it =  tmp.begin() +(rand() % ((int) tmp.size()) );
         spaces.push_back( *(it) );
         tmp.erase(it);

    }

    return spaces;
}


//[[Rcpp::export]]

vector<int> fromWhichHowMany(int minContigLength,int totalLength,vector<int>& lengths,int needed){
    default_random_engine generator;
    vector<int> res;
    res.reserve(lengths.size());
    int share;
    int tmp;
    for(vector<int>::iterator it = lengths.begin();it != lengths.end();it++){
        if((*it) >= minContigLength){
            share = round(needed *((*it)/(double)totalLength));
            poisson_distribution<int> distribution(share);
            tmp = distribution(generator) %(*it);
            distribution.reset();
            needed -= tmp;
            totalLength -= (*it);
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

list<list<list<vector<int> > > > mkContigs(list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int seed = 0,bool sameContigLength = true){

    srand(seed);
    list<list<list<vector<int> > > > res;
    int at = 0;
    vector<int>::iterator totLen = lengthSums.begin();
    
    for(list<vector<int> >::iterator it = lengths.begin(); it != lengths.end();it++){
        
        list<list<vector<int> > > tmp;
        vector<int> baseNrs;
        vector<int> indicies;
        int count = 0;
        
        double partCovered = (1+ (rand() % 9))/10.0;
        baseNrs.push_back((int) ((*totLen) *partCovered));
        double compPart = (6+ (rand() % 4))/10.0;
        baseNrs.push_back((*prev(baseNrs.end()) / (compPart * 100.0)*100 -(*prev(baseNrs.end()))));
        
        vector<int> bigEnough;
        for(int n = 0;n < (int) lengthSums.size();n++){
            if(n != count && lengthSums[n] > (*prev(baseNrs.end()))){
                bigEnough.push_back(n);
            }
        }
        
        int index = bigEnough[rand() % (bigEnough.size() -1)];
        indicies.push_back(count);
        indicies.push_back(index);
        vector<int> chromBaseNrs;
        
        for(int n = 0; n < (int) baseNrs.size();n++){
            
            list<vector<int> > tmp2;
            chromBaseNrs = fromWhichHowMany(minContigLength,(*totLen),(*next(it,indicies[n])),baseNrs[n]);
            tmp.push_back(list<vector<int> > {vector<int> (1,indicies[n])});
            for(int j = 0; j < (int) (*it).size();j++){
                
                vector<int> starts;
                vector<int> ends;
                vector<int> contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[j]);
                vector<int> spaces = randomSpaces((int)contigs.size() +1,(*it)[j] -chromBaseNrs[j]);
                vector<int>::iterator co = contigs.begin();
                vector<int>::iterator sp = spaces.begin();
                
                for(int n = 0;n < (int)contigs.size();n++){
                    at += (*sp);
                    starts.push_back(at+1);
                    at += (*co);
                    ends.push_back(at);
                    sp++;
                    co++;
                }
                
                at = 0;
                tmp2.push_back(starts);
                tmp2.push_back(ends);
            }
            tmp.push_back(tmp2);
            count++;
        }
        res.push_back(tmp);
    }
    return res;
}
