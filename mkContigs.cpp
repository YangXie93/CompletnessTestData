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
#include <string>

using namespace std;

//[[Rcpp::export]]

vector<int> randomContigs(int minContigLength,int meanContigLength,int covering,string distr = "normal",int seed = 0){
    if(covering <= minContigLength){
        return vector<int> (1,covering);
    }
    default_random_engine generator;
    generator.seed(seed);
    normal_distribution<double> normDist(meanContigLength,sqrt(meanContigLength));
    poisson_distribution<int> poisDist(meanContigLength);
    exponential_distribution<double> expDist(meanContigLength);
    uniform_int_distribution<int> uniformDist(minContigLength,meanContigLength*2);
    vector<int> res;
    int tmp = 0;
    while(covering > 0){
        if(distr == "poisson" || distr == "normal" || distr == "exponential" || distr == "uniform"){
            if(distr == "poisson"){
                tmp = (int) poisDist(generator);
                poisDist.reset();
            }
            if(distr == "normal"){
                tmp = (int) normDist(generator);
                normDist.reset();
            }
            if(distr == "exponential"){
                tmp = (int) expDist(generator);
                expDist.reset();
            }
            if(distr == "uniform"){
                tmp = (int) uniformDist(generator);
                uniformDist.reset();
            }
        }
        else{
            Rcpp::Rcerr << "die gewählte verteilung steht nicht zur wahl." << endl << "verfügbar sind: \"normal\", \"poisson\",\"exponential\" und \"uniform\" " << endl;
            Rcpp::Rcerr << "Das Programm wird jetzt mit einer Normal verteilung fortgesetzt." << endl;
            tmp = (int) normDist(generator);
            normDist.reset();
        }
        if(tmp >= covering){
            res.push_back(covering);
            covering = 0;
        }
        else{
            if(tmp > 0){
                covering -= tmp;
                res.push_back(tmp);
            }
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

vector<int> fromWhichHowMany(int minContigLength,int totalLength,vector<int> lengths,int needed){
    default_random_engine generator;
    vector<int> res;
    res.reserve(lengths.size()/2);
    int share;
    int tmp;
    vector<int> tmp1;
    vector<int> accession;
    tmp1.reserve(lengths.size()/2);
    accession.reserve(lengths.size()/2);
    int j = 0;
    for(int i = 1;i < (int) lengths.size();i+= 2){
        res.push_back(0);
        tmp1.push_back(lengths[i]);
        accession.push_back(j);
        j++;
    }
    for(int i = 0; i < (int) res.size();i++){
        if(needed > 0){
            int x = (rand() % (tmp1.size()));
            vector<int>::iterator it = next(tmp1.begin(),x);
            vector<int>::iterator acc = next(accession.begin(),x);
            share = round(needed *((*it)/(double)totalLength));
            poisson_distribution<int> distribution(share);
            tmp = distribution(generator);
            if(tmp > (*it)){
                tmp -= tmp %(*it); 
            }
            if(tmp >= minContigLength){
                needed -= tmp;
                totalLength -= (*it);
                res[(*acc)] = tmp;
            }
            else{
                needed -= minContigLength;
                totalLength -= (*it);
                res[(*acc)] = minContigLength;
            }
            tmp1.erase(it);
            accession.erase(acc);
        }
    }

    return res;
}


//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::export]]

list<list<list<vector<int> > > > mkContigs(list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,int seed = 0,string distr = "normal"){

    default_random_engine generator;
    generator.seed(seed);
    srand(seed);
    list<list<list<vector<int> > > > res;
    int at = 0;

    
    for(int i = 0; i < number; i++){
        vector<int> bigEnough;
        double partCovered = (60+ (rand() % 40))/100.0;
        for(int n = 0;n < (int) lengthSums.size();n++){
            if((lengthSums[n] *partCovered) >= minContigLength){
                bigEnough.push_back(n);
            }
        }
        int which = generator() % bigEnough.size();
        vector<int>::iterator totLen = next(lengthSums.begin(),bigEnough[which]);
        list<list<vector<int> > > tmp;
        vector<int> baseNrs;
        vector<int> indicies;


        baseNrs.push_back((int) ((*totLen) *partCovered));
        double contPart = ((rand() % 40))/100.0;
        if((*totLen) *contPart <= minContigLength){
            baseNrs.push_back(minContigLength);
        }
        else{
            baseNrs.push_back((*totLen) *contPart);
        }
        bigEnough.clear();
        for(int n = 0;n < (int) lengthSums.size();n++){
            if(n != which && lengthSums[n] >= (*prev(baseNrs.end()))){
                bigEnough.push_back(n);
            }
        }
        if((int) bigEnough.size() != 0){
            int index = bigEnough[(rand() % bigEnough.size()) -1];
            indicies.push_back(distance(lengthSums.begin(),totLen));
            indicies.push_back(index);
            vector<int> chromBaseNrs;
            for(int n = 0; n < (int) baseNrs.size();n++){
    
                list<vector<int> > tmp2;
                chromBaseNrs = fromWhichHowMany(minContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n]);
                int l = 0;
                for(int j = 1; j < (int) (*next(lengths.begin(),indicies[n])).size();j += 2){
                    vector<int> contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[l],distr,seed);
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
        else{
            i--;
         }
    }
    return res;
}

