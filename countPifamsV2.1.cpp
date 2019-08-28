#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <Rcpp.h>
#include <time.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <string>
#include <cmath>

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
    geometric_distribution<int> expDist;
    uniform_int_distribution<int> uniformDist(minContigLength,meanContigLength*2);
    vector<int> res;
    int tmp = 0;
    int mean = 0;
    int isSmalerZero = 0;
    while(covering > 0){
        if(distr == "poisson" || distr == "normal" || distr == "exponential" || distr == "uniform"){
            if(distr == "poisson"){
                tmp = poisDist(generator);
                poisDist.reset();
            }
            if(distr == "normal"){
                tmp = (int) normDist(generator);
                normDist.reset();
            }
            if(distr == "exponential"){

                tmp = minContigLength + (expDist(generator));
                if(tmp > minContigLength){
                    isSmalerZero++;
                }
                expDist.reset();
            }
            if(distr == "uniform"){
                tmp = uniformDist(generator);
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
            covering = 0;
        }
        else{
            if(tmp > 0){
                covering -= tmp;
                res.push_back(tmp);
                mean += tmp;
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

vector<vector<int> > mkContigs(list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,vector<int>& pseuPoint,vector<double> comp,vector<double> cont,int seed = 0,string distr = "normal"){

    default_random_engine generator;
    generator.seed(seed);
    srand(seed);
    
    vector<vector<int> > res;
    res.reserve(number);

    vector<int>::iterator totLen;
    vector<int> baseNrs;
    vector<int> indicies;
    vector<int> contigs;
    vector<int> spaces;
    vector<int>::iterator co;
    vector<int>::iterator sp;
    vector<int> starts;
    vector<int> ends;
    vector<int> chromBaseNrs;
    
    int at = 0;
    double partCovered;
    vector<int>::iterator which;
    bool swtch;
    int index;
    int l;
    int n;
    int j;
    int m;
    int count;
    
    if(cont[1] > comp[0]){
        cont[1] = comp[0];
    }
    
    
    for(int i = 0; i < number; i++){
        partCovered = ((comp[0] *100) + (rand() % (int)((comp[1]-comp[0]) *100)))/100.0;
        which = next(lengthSums.begin(),(generator() % lengthSums.size()));
        
        count = 0;
        while(((*which)* partCovered) < minContigLength && count < lengthSums.size()){
            which++;
            count++;
            if(which == lengthSums.end()){
                which = lengthSums.begin();
            }
        }
        if(count == lengthSums.size() && (*which) < minContigLength){
            Rcpp::Rcerr << "Die mindest Länge ist zu groß für den Datensatz" << endl;
            return res;
        }
        
        totLen = which;

        baseNrs.push_back((int) ((*totLen) *partCovered));
        double contPart = ((cont[0] *100) + (rand() % (int)(cont[1] *100)))/100.0;
        if((*totLen) *contPart <= minContigLength){
            baseNrs.push_back(minContigLength);
        }
        else{
            baseNrs.push_back((*totLen) *contPart);
        }
        
        count = 0;
        which = next(lengthSums.begin(),(generator() % lengthSums.size()));
        while((*which) < *prev(baseNrs.end()) && count < lengthSums.size()){
            which++;
            count++;
            if(which == lengthSums.end()){
                which = lengthSums.begin();
            }
            if(which == totLen){
                which++;
            }
        }
        if((*which) < *prev(baseNrs.end()) && count == lengthSums.size()){
            Rcpp::Rcerr << "Die geforderte Completeness ist zu hoch für den Datensatz" << endl;
            return res;
        }
        res.push_back(vector<int> (0));
        index = distance(lengthSums.begin(),totLen);
        indicies.push_back(distance(lengthSums.begin(),totLen));
        indicies.push_back(index);
        pseuPoint.push_back(distance(lengthSums.begin(),totLen));
        pseuPoint.push_back(index);
        for(n = 0; n < (int) baseNrs.size();n++){

            res.push_back(vector<int> (0));
            chromBaseNrs = fromWhichHowMany(minContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n]);
            l = 0;
            for(j = 1; j < (int) (*next(lengths.begin(),indicies[n])).size();j += 2){
                contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[l],distr,seed);
                spaces = randomSpaces((int)contigs.size() +1,(*next(lengths.begin(),indicies[n]))[j] -chromBaseNrs[l]);
                co = contigs.begin();
                sp = spaces.begin();

                starts.reserve(contigs.size());

                ends.reserve(contigs.size());
                swtch = true;
                for(int m = 0;m < (int)contigs.size();m++){
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
                res.push_back(vector<int> {(*next(lengths.begin(),indicies[n]))[j-1]});
                res.push_back(starts);
                res.push_back(ends);
                starts.clear();
                ends.clear();
                l++;
            }
        }
        baseNrs.clear();
        indicies.clear();
    }
    return res;
}


//[[Rcpp::export]]

vector<int> intervallOverlap(int start1, int end1,int start2,int end2){
    int s;
    int e;
    if(start2 > end1 || end2 < start1 ){
        return {0};
    }
    if(start1 > start2){
        s = start1;
        if(end2 < end1){
            e = end2;
        }
        else{
            e = end1;
        }
    }
    else{
        s = start2;
        if(end2 < end1){
            e = end2;
        }
        else{
            e = end1;
        }
    }
    return {s,e};
}

//[[Rcpp::export]]

Rcpp::List countPifams(list<list<vector<int> > >& pifams,list<list<vector<int> > >& ORFs,vector<vector<int> > contigs,vector<int>& names){

    Rcpp::List res;
    vector<int> id;
    vector<int> times;
    vector<int> baseNum;
    vector<bool> notInOrf;
    vector<int> name;
    vector<double> completeness;
    vector<double> contamination;

    int o = 0;
    int n = 0;
    int x = 0;
    int zwsch = 0;
    int total = 0;
    int partial = 0;
    int contTotal = 0;
    int transfer;
    int GenNr;

    vector<int> tmp1;

    vector<int>::iterator trans;
    vector<int>::iterator t;
    vector<int>::iterator starts;
    vector<int>::iterator ends;
    vector<int>::iterator width;
    vector<int>::iterator values;
    vector<int>::iterator orfValues;
    vector<int>::iterator orfWidths;

    vector<int>::iterator nms = names.begin();
    list<vector<int> >::iterator orfs;
    list<vector<int> >::iterator pis;

    vector<vector<int> >::iterator con = contigs.begin();

    vector<int>::iterator bN;
    vector<bool>::iterator nio;
    vector<int>::iterator tms;
    bool swtch = false;

    while(distance(con,contigs.end()) > 0){
        // alle Werte
        Rcpp::List tmp;
        swtch = false;
        GenNr = 1;
        while((int)(*con).size() != 0 && distance(con,contigs.end()) > 0){
            swtch = true;
            // alle Genome
            while((int) (*con).size() != 0 && distance(con,contigs.end()) > 0){
                // alle Chromosomen

                name.push_back(*((*con).begin()));
                transfer = (*nms);
                orfs = (*next(ORFs.begin(),transfer)).begin();
                for(pis = (*next(pifams.begin(),transfer)).begin();pis != (*next(pifams.begin(),transfer)).end();pis = next(pis,2)){
                    // für alle 6 GENOME und ORF Werte

                    starts = (*next(con)).begin();
                    ends = (*next(con,2)).begin();
                    width = (*pis).begin();
                    values = (*next(pis)).begin();
                    orfWidths = (*orfs).begin();
                    orfValues = (*next(orfs)).begin();

                    for(starts = (*next(con)).begin();starts != (*next(con)).end();starts++){
                        //für alle Length Werte in ORF
                        x += (*ends) - (*starts) +1;
                        while((*starts) > (o +(*orfWidths))){
                            o += (*orfWidths);
                            orfWidths++;
                            orfValues++;
                        }
                        while(o > (*starts) || (o + (*orfWidths)) < (*ends)){
                            if((*orfValues) > 0){
                                fill(notInOrf.begin(),notInOrf.end(),true);
                                while(o > n+(*width)){
                                    n += (*width);
                                    width++;
                                    values++;
                                }
                                while(n > (*starts) || (n + (*width)) < (*ends)){
                                    if((*values) > 0){
                                        tmp1 = intervallOverlap(n+1,n+(*width),(*starts),(*ends));
                                        if((int) tmp1.size() > 1){
                                            t = find(id.begin(),id.end(),(*values));
                                            if(t == id.end()){
                                                baseNum.push_back(tmp1[1] - tmp1[0] +1);
                                                times.push_back(1);
                                                id.push_back((*values));
                                                notInOrf.push_back(false);
                                            }
                                            else{
                                                bN = next(baseNum.begin(),t-id.begin());
                                                nio = next(notInOrf.begin(),t-id.begin());
                                                if((*nio)){
                                                    tms = next(times.begin(),t-id.begin());
                                                    (*tms)++;
                                                    (*nio) = false;
                                                }
                                                (*bN) += tmp1[1] - tmp1[0] +1;
                                            }
                                        }
                                    }
                                    if((n + (*width)) < (*ends)){
                                        n += (*width);
                                        width++;
                                        values++;
                                    }
                                    else{
                                        break;
                                    }
                                }
                            }
                            if((o + (*orfWidths)) < (*ends)){
                                o += (*orfWidths);
                                orfWidths++;
                                orfValues++;
                            }
                            else{
                                break;
                            }
                        }

                        ends++;
                    }
                    //end alle Length Werte aus ORF
                    if(pis == (*next(pifams.begin(),transfer)).begin()){
                        total += o;
                        while(orfWidths != (*orfs).end()){
                            total += (*orfWidths);
                            orfWidths++;
                        }
                    }
                    if(x > zwsch){
                        zwsch = x;
                    }
                    o = 0;
                    n = 0;
                    x = 0;
                    orfs = next(orfs,2);
                }
                //ende alle 6 GENOME und ORF Werte
                if(GenNr == 1){
                    partial += zwsch;
                }
                else{
                    contTotal += zwsch;
                }
                zwsch = 0;
                con = next(con,3);
                nms++;
            }
            //ende alle Chromosomen
            if(GenNr == 1){
                completeness = {partial/(double) total};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("compChromID") = name,Rcpp::Named("compPifamNames") = id,Rcpp::Named("compPifamCount") = times,Rcpp::Named("compBaseCount") = baseNum,Rcpp::Named("completness") = completeness));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
            else{
                contamination = {(contTotal/(double)total)};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("contChromID") = name,Rcpp::Named("contPifamNames") = id,Rcpp::Named("contPifamCount") = times,Rcpp::Named("contBaseCount") = baseNum,Rcpp::Named("contamination") = contamination));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
            total = 0;
            partial = 0;
            contTotal = 0;
            con++;
            GenNr++;
        }
        //ende alle Genome
        if(swtch){
            res.push_back(tmp);
        }
        con++;
    }
    return res;
}


//[[Rcpp::export]]

Rcpp::List compTestData(list<list<vector<int> > >& pifams,list<list<vector<int> > >& ORFs,list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,vector<double> comp,vector<double> cont,int seed = 0,string distr = "normal"){
    vector<int>names;

    return countPifams(pifams,ORFs,mkContigs(lengths,lengthSums,minContigLength,meanContigLength,number,names,comp,cont,seed,distr),names);
}
