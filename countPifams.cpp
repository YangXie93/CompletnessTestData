#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <Rcpp.h>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <algorithm>
#include <string>

using namespace std;

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
    Rcpp::Rcout << "RandCont: how many: " << res.size() << endl;
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
    for(vector<int>::iterator it = next(lengths.begin());distance(it,lengths.end()) > 0;it = next(it,2)){
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
    Rcpp::Rcout << "FWHM: needed: " << needed << " Nr.: " << res[0] << endl;
    
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
        Rcpp::Rcout << "At: " << i << endl;
        int which = generator() % lengths.size();
        vector<int>::iterator totLen = next(lengthSums.begin(),which);
        list<list<vector<int> > > tmp;
        vector<int> baseNrs;
        vector<int> indicies;
        
        double partCovered = (60+ (rand() % 40))/100.0;
        Rcpp::Rcout << "completeness: " << partCovered << endl;
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
        Rcpp::Rcout << "contamination: " << 1- compPart << endl;
        if((*prev(baseNrs.end()) / (compPart * 100.0)*100 -(*prev(baseNrs.end()))) <= minContigLength){
            Rcpp::Rcout << "entschuldigt!" << endl;
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
        Rcpp::Rcout<< "index1: "  << which << " index2: " << index << endl;
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
        Rcpp::Rcout << endl;
        res.push_back(tmp);
    }
    Rcpp::Rcout << "mkContigs end" << endl;
    return res;
}



//[[Rcpp::export]]

vector<int> intervallOverlap(int start1, int end1,int start2,int end2){
    vector<int> res;
    if(start2 > end1 || end2 < start1 ){
        res.push_back(0);
        return res;
    }
    if(start1 > start2){
        res.push_back(start1);
        if(end2 < end1){
            res.push_back(end2);
        }
        else{
            res.push_back(end1);
        }
    }
    else{
        res.push_back(start2);
        if(end2 < end1){
            res.push_back(end2);
        }
        else{
            res.push_back(end1);
        }
    }
    return res;
}

//[[Rcpp::export]]

Rcpp::List countPifams(list<list<vector<int> > >& pifams,list<list<vector<int> > >& ORFs,list<list<list<vector<int> > > >& contigs,vector<int>& names){

    Rcpp::List res;
    vector<int> id;
    vector<int> times;
    vector<int> baseNum;
    vector<bool> notInOrf;
    vector<int> name;
    
    for(list<list<list<vector<int> > > >::iterator con = contigs.begin();con != contigs.end();con++){
       
        Rcpp::List tmp;
        int total = 0; 
        int partial = 0;
        int contTotal = 0;
        
        for(list<list<vector<int> > >::iterator cons = (*con).begin();cons != (*con).end();cons++){
            
            for(list<vector<int> >::iterator conts = (*cons).begin();conts != (*cons).end();conts = next(conts,3)){
                int o = 0;
                int n = 0;
                int x = 0;
                name.push_back(*((*conts).begin()));
                vector<int>::iterator trans = find(names.begin(),names.end(),*((*conts).begin()));
                if(trans != names.end()){
                    int transfer = distance(names.begin(),trans);
                    
                    for(vector<int>::iterator wi = (*(*next(pifams.begin(),transfer)).begin()).begin();wi != (*(*next(pifams.begin(),transfer)).begin()).end();wi++){
                        total += (*wi);
                    }
                    vector<int>::iterator e = (*next(conts,2)).begin();
                    int t = 0;
                    for(vector<int>::iterator s = (*next(conts)).begin();s != (*next(conts)).end();s++){
                        t += (*e)-(*s) +1;
                        e++;
                    }
                    if(cons == (*con).begin()){
                        partial += t;
                    }
                    contTotal += t;
                    list<vector<int> >::iterator orfs = (*next(ORFs.begin(),transfer)).begin();
                    
                    for(list<vector<int> >::iterator pis = (*next(pifams.begin(),transfer)).begin();pis != (*next(pifams.begin(),transfer)).end();pis = next(pis,2)){
                        // für alle 6 GENOME und ORF Werte
                        
                        vector<int>::iterator starts = (*next(conts)).begin();
                        vector<int>::iterator ends = (*next(conts,2)).begin();
                        vector<int>::iterator width = (*pis).begin();
                        vector<int>::iterator values = (*next(pis)).begin();
                        vector<int>::iterator orfValues = (*next(orfs)).begin();
                        
                        bool swtch = true;
                        for(vector<int>::iterator orfWidths = (*orfs).begin(); distance(orfWidths,(*orfs).end()) > 0; orfWidths++){
                            //für alle Length Werte in ORF
                            if(starts != (*next(conts)).end()){
                                for(int i = 0;i < (int) notInOrf.size();i++){
                                    notInOrf[i] = true;
                                }
                                while(width != (*pis).end() && o+(*orfWidths) > n && swtch){
                                
                                    
                                    if((*orfValues) >= 0){
                                        while(starts != (*next(conts)).end() && n+(*width) > (*starts) && swtch){
                                            
                                            if((*values) >= 0){
                                                vector<int> tmp = intervallOverlap(n,n+(*width),(*starts),(*ends));
                                                if((int) tmp.size() > 1){
                                                    vector<int>::iterator t = find(id.begin(),id.end(),(*values));
                                                    if(t == id.end()){
                                                        baseNum.push_back(tmp[1] - tmp[0] +1);
                                                        times.push_back(1);
                                                        id.push_back((*values));
                                                        notInOrf.push_back(false);
                                                    }
                                                    else{
                                                        if(notInOrf[((int) (t-id.begin()))]){
                                                            times[((int) (t-id.begin()))]++;
                                                            notInOrf[((int) (t-id.begin()))] = false;
                                                        }
                                                        baseNum[((int) (t-id.begin()))] += tmp[1] - tmp[0];
                                                    }
                                                }
                                            }
                                            if((*ends) >= (n+(*width))){
                                                swtch = false;
                                            }
                                            else{
                                                x += (*ends) -(*starts);
                                                starts++;
                                                ends++;
                                            }
                                        }
                                        
                                        swtch = true;
                                    }

                                    n += (*width);
                                    values++;
                                    width++;
                                    
                                }
                                swtch = true;
                            }
                            o += (*orfWidths);
                            orfValues++;
                        }
                        //end alle Length Werte aus ORF
                        o = 0;
                        n = 0;
                        x = 0;
                        orfs = next(orfs,2);
                    }
                    //ende alle 6 GENOME und ORF Werte
                    
                }
                if(cons == (*con).begin()){
                    tmp.push_back(Rcpp::List::create(Rcpp::Named("compChromID") = name,Rcpp::Named("compPifamNames") = id,Rcpp::Named("compPifamCount") = times,Rcpp::Named("compBaseCount") = baseNum));
                    name.clear();
                    id.clear();
                    baseNum.clear();
                    times.clear();
                    notInOrf.clear();
                }
                else{
                    tmp.push_back(Rcpp::List::create(Rcpp::Named("contChromID") = name,Rcpp::Named("contPifamNames") = id,Rcpp::Named("contPifamCount") = times,Rcpp::Named("contBaseCount") = baseNum));
                    name.clear();
                    id.clear();
                    baseNum.clear();
                    times.clear();
                    notInOrf.clear();
                }
            }
            //ende alle Chromosomen

            if(cons == (*con).begin()){
                vector<double> completeness = {(double)partial/(double) total};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("completness") = completeness));
            }
            else{
                vector<double> contamination = {1 -((double)partial/(double)contTotal)};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("contamination") = contamination));
            }
        
        }
        
        //ende alle Genome
        res.push_back(tmp);
    }
    Rcpp::Rcout << "countPifams end" << endl;
    return res;
}

//[[Rcpp::export]]

Rcpp::List mkContigsAndCountPifams(list<vector<int> >& lengths,vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,list<list<vector<int> > >& pifams,list<list<vector<int> > >& ORFs,vector<int>& names,int seed = 0,string distr = "normal"){
    list<list<list<vector<int> > > > contigs = mkContigs(lengths,lengthSums,minContigLength,meanContigLength,number,seed,distr);
    Rcpp::List res = countPifams(pifams,ORFs,contigs,names);
    return res;
}


