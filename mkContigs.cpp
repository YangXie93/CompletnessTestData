#include <iostream>
#include <vector>
#include <list>
#include <Rcpp.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <algorithm>

using namespace std;

//[[Rcpp::export]]

vector<int> randomContigs(int minContigLength,int meanContigLength,int covering){
    if(covering <= minContigLength){
        return vector<int> (1,covering);
    }
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
    if(free <= 0){
        return vector<int> (numSp,0);
    }
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
    for(vector<int>::iterator it = next(lengths.begin());it != lengths.end();it++){
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

list<list<list<vector<int> > > > mkContigs(list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,int seed = 0){

    srand(seed);
    list<list<list<vector<int> > > > res;
    int at = 0;

    for(int i = 0; i < number; i++){
        int which = rand() % (int) lengths.size();
        vector<int>::iterator totLen = next(lengthSums.begin(),which);
        list<list<vector<int> > > tmp;
        vector<int> baseNrs;
        vector<int> indicies;

        double partCovered = (60+ (rand() % 40))/100.0;
        if((*totLen)* partCovered < minContigLength){
            if((*totLen) < minContigLength){
                baseNrs.push_back((*totLen));
            }
            else{
                baseNrs.push_back(minContigLength);
            }
        }
        else{
            baseNrs.push_back((int) ((*totLen) *partCovered));
        }
        double compPart = (60+ (rand() % 40))/100.0;
        if((*prev(baseNrs.end()) / (compPart * 100.0)*100 -(*prev(baseNrs.end()))) <= minContigLength){
            baseNrs.push_back(minContigLength);
        }
        else{
            baseNrs.push_back((*prev(baseNrs.end()) / (compPart * 100.0)*100 -(*prev(baseNrs.end()))));
        }
        vector<int> bigEnough;
        for(int n = 0;n < (int) lengthSums.size();n++){
            if(n != which && lengthSums[n] > (*prev(baseNrs.end()))){
                bigEnough.push_back(n);
            }
        }

        int index = bigEnough[rand() % (bigEnough.size() -1)];
        indicies.push_back(which);
        indicies.push_back(index);
        vector<int> chromBaseNrs;
        for(int n = 0; n < (int) baseNrs.size();n++){

            list<vector<int> > tmp2;
            chromBaseNrs = fromWhichHowMany(minContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n]);
            int l = 0;
            for(int j = 1; j < (int) (*next(lengths.begin(),indicies[n])).size();j += 2){
                vector<int> contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[l]);
                vector<int> spaces = randomSpaces((int)contigs.size() +1,(*next(lengths.begin(),indicies[n]))[j] -chromBaseNrs[l]);
                vector<int>::iterator co = contigs.begin();
                vector<int>::iterator sp = spaces.begin();
                vector<int> starts;
                starts.reserve(contigs.size());
                vector<int> ends;
                ends.reserve(contigs.size());
                bool swtch = true;
                for(int n = 0;n < (int)contigs.size();n++){
                    at += (*sp);
                    if(swtch){
                        starts.push_back(at+1);
                    }
                    at += (*co);
                    sp++;
                    co++;
                    if((*sp) > 0){
                        ends.push_back(at);
                    }
                    else{
                      swtch = false;
                    }
                }
                if(ends.size() == 0){
                    ends.push_back(at);
                }

                at = 0;
                tmp2.push_back(vector<int> {(*next(lengths.begin(),indicies[n]))[j-1]});
                tmp2.push_back(starts);
                tmp2.push_back(ends);
                l++;
            }
            tmp.push_back(tmp2);
        }
        res.push_back(tmp);
    }
    return res;
}
