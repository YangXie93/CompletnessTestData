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

//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;



std::vector<int> randomContigs(int minContigLength,int meanContigLength,int covering,std::string distr = "normal",int seed = 0){
    if(covering <= 0){
        return std::vector<int> {};
    }
    if(covering <= minContigLength){
        return std::vector<int> (1,covering);
    }
    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> normDist(meanContigLength,sqrt(meanContigLength));
    std::poisson_distribution<int> poisDist(meanContigLength);
    std::geometric_distribution<int> expDist;
    std::uniform_int_distribution<int> uniformDist(minContigLength,meanContigLength*2);
    std::vector<int> res;
    int tmp = 0;
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
                expDist.reset();
            }
            if(distr == "uniform"){
                tmp = uniformDist(generator);
                uniformDist.reset();
            }
        }
        else{
            Rcerr << "die gewählte verteilung steht nicht zur wahl." << std::endl << "verfügbar sind: \"normal\", \"poisson\",\"exponential\" und \"uniform\" " << std::endl;
            Rcerr << "Das Programm wird jetzt mit einer Normal verteilung fortgesetzt." << std::endl;
            tmp = (int) normDist(generator);
            normDist.reset();
        }
        if(tmp >= covering){
            if(covering >= minContigLength){
                res.push_back(covering);
            }
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

//' @export
//[[Rcpp::export]]
std::vector<int> randomSpaces(int numSp,int free,int seed){
    if(free <= 0){
        return std::vector<int> (numSp,0);
    }
    std::vector<int > spaces;
    spaces.reserve(numSp);

    std::default_random_engine generator;
    generator.seed(seed);

    std::vector<int> t;
    t.reserve(numSp);
    double sum = 0.0;
    int flex;
    int i;
    int sum2 = 0;
    int* elem = t.data();

    for(i = 0;i < numSp;i++){
        flex = generator() % 10000;
        t.push_back(flex);
        sum += flex;
    }

    for(i = 0;i < numSp;i++){
        *(elem +i) =  (int)(free*(*(elem +i)/sum));
        sum2 += *(elem +i);
    }
    if(sum2 > free){
        while(sum2 > free){
            i = generator() % numSp;
            if(*(elem +i) > 0){
                (*(elem +i))--;
                sum2--;
            }
            i++;
        }
    }
    if(sum2 < free){
        while(sum2 < free){
            i = generator() % numSp;
            (*(elem +i))++;
            sum2++;
        }
    }
    std::vector<int>::iterator it;
    while((int) t.size() > 0){
         it =  t.begin() +(generator() % ((int) t.size()));
         spaces.push_back( *(it));
         t.erase(it);
    }
    return spaces;
}

int whichToSmall;

//[[Rcpp::export]]
std::vector<int> fromWhichHowMany(int minContigLength,int totalLength,std::vector<int> lengths,int needed,int seed){
    
    std::default_random_engine generator;
    generator.seed(seed);
    std::vector<int> res (lengths.size(),0);
    int share;
    int at;
    int tmp;
    int ges;
    int min;
    int relMin;
    
    std::vector<int> tmp1;
    std::vector<int> accession;
    std::vector<int> acc;
    
    int j = 0;
    int n = 0;
    for(int i = 0;i < (int) lengths.size();i++){
        if(lengths[i] > minContigLength){
            tmp1.push_back(lengths[i]);
            accession.push_back(j);
            acc.push_back(n);
            n++;
        }
        else{
            totalLength -= lengths[i];
        }
        j++;
    }
    
    if(tmp1.size() > 0){
        
        std::vector<int> tmp2;
        
        for(int i = 0;i < tmp1.size();i++){
            relMin = (int) round(((needed/(double)totalLength) * (tmp1[i]/(double)totalLength)) *totalLength);
            min = (int) round(((minContigLength)/(double) (needed - minContigLength)) *(totalLength-tmp1[i]));
            
            if(min < tmp1[i]){
                if(relMin > min){
                    tmp2.push_back(relMin + (generator() % (tmp1[i] - relMin) ));
                }
                else{
                    tmp2.push_back(min + (generator() % (tmp1[i] - min)));
                }
            }
            else{
                tmp2.push_back(0);
            }
        }
        ges = accumulate(tmp2.begin(),tmp2.end(),0);
        
        if(ges == 0){
            int x = 0;
            int max = 0;
            int nenn = 1;
            while(tmp1[x] < needed){
                x++;
                if(tmp1[x] > tmp1[max]){
                    max = x;
                }
                if(x == tmp1.size() -1){
                    x = max;
                    nenn = tmp1[x];
                    break;
                }
            }
            tmp2[x] = nenn;
            ges = needed;
        }
        
        for(int i = 0;i < tmp2.size();i++){
            res[accession[i]] = (int) (needed* (tmp2[i]/(double)ges));
            // if(res[accession[i]] < minContigLength){
            //     Rcout << whichToSmall << ": kleiner als minContigLength \n";
            //     Rcout << needed << " " << tmp2[i]<< " " << ges << " " << res[accession[i]] << " " << tmp1[i] << " " << tmp1.size() << std::endl;
            // }
        }
        
    }
    return res;
}

int minL;

bool isBiggerMinL(int i){
    return(i > minL);
}


//' @export
//[[Rcpp::export]]
std::list<std::list<std::list<std::vector<int> > > > mkContigs(std::list<std::vector<int> >& lengths,std::list<std::vector<int> > &IDs,std::vector<int>& lengthSums,int minContigLength,int meanContigLength,int number,std::vector<double>& comp,std::vector<double>& cont,int seed = 0,std::string distr = "normal"){
    whichToSmall = 0;
    minL = minContigLength;
    
    std::default_random_engine generator;
    generator.seed(seed);
    
    std::list<std::list<std::list<std::vector<int> > > > res;
    std::vector<int>::iterator totLen;
    std::vector<int> baseNrs;
    std::vector<int> indicies;
    std::vector<int> contigs;
    std::vector<int> spaces;
    std::vector<int>::iterator co;
    std::vector<int>::iterator sp;
    std::vector<int> starts;
    std::vector<int> ends;
    std::vector<int> chromBaseNrs;
    std::vector<int>::iterator which;
    std::vector<int>::iterator max;
    std::list<std::vector<int> >::iterator tester;
    
    bool swtch;
    int index;
    int l;
    int n;
    int j;
    int count;
    int accuContigs;
    int at = 0;
    int contigSum;
    double partCovered;
    bool contIsNull;
    bool justZero = false;
    int testTmp;
    
    //--------- Prüfen von comp und cont ---------------
    
    if(cont[1] > comp[0]){
        cont[1] = comp[0];
    }
    if(cont[0] == cont[1]){
        cont[0] = cont[1] -0.01;
    }
    if(comp[0] == comp[1]){
        comp[0] = comp[1] -0.01;
    }
    
    //--------------------------------------------------
    
    for(int i = 0; i < number; i++){
        
        whichToSmall++;
        count = 0;
        contIsNull = false;
        
        //------------------------ Aussuchen des bin angebenden Genoms und der completness ----------------------
        
        partCovered = ((comp[0] *100) + (generator() % (int)((comp[1]-comp[0]) *100)))/100.0;
        testTmp = (generator() % lengthSums.size());
        which = next(lengthSums.begin(),testTmp);
        tester = next(lengths.begin(),testTmp);
        
        while(((*which)* partCovered) < minContigLength && count < (int)lengthSums.size()){
            which++;
            count++;
            tester++;

            if((find_if((*tester).begin(),(*tester).end(),isBiggerMinL)) == (*tester).end()){
                which++;
                tester++;
                count++;
            }
            
            if(which == lengthSums.end()){
                which = lengthSums.begin();
            }
        }
        
        if(count == (int)lengthSums.size() && (*which) < minContigLength){
            Rcerr << "Die mindest Länge ist zu groß für den Datensatz" << std::endl;
            return res;
        }
        totLen = which;
        baseNrs.push_back((int) ((*totLen) *partCovered));
        
        //-------------------------------------------------------------------------------------------------------
        
        //------------------------- Aussuchen des genoms zur contamination --------------------------------------
        
        double contPart = ((cont[0] *100) + (generator() % ((int)(cont[1]*100) - (int)(cont[0] *100))))/100.0;
        if((*totLen) *contPart <= minContigLength){
            if(contPart != 0){
                baseNrs.push_back(minContigLength);
            }
            else{
                contIsNull = true;
            }
        }
        else{
            baseNrs.push_back((*totLen) *contPart);
        }
        
        if(!contIsNull){
            count = 0;
            testTmp = (generator() % lengthSums.size());
            which = next(lengthSums.begin(),testTmp);
            tester = next(lengths.begin(),testTmp);
            max = lengthSums.begin();
            while(((*which) < *prev(baseNrs.end()) && count < (int)lengthSums.size()) || which == totLen){
                if(which == lengthSums.end()){
                    which = lengthSums.begin();
                }
                else{
                    
                    if(which > max && which != totLen){
                        max = which;
                    }
                    which++;
                    count++;
                    tester++;

                    if((find_if((*tester).begin(),(*tester).end(),isBiggerMinL)) == (*tester).end()){
                        which++;
                        tester++;
                        count++;
                    }
                }
            }
            if(((*which) < *prev(baseNrs.end()) && count == (int)lengthSums.size()) || which == lengthSums.end() || which == totLen){
                baseNrs.pop_back();
            }
        }
        
        //--------------------------------------------------------------------------------------------------------
        
        index = distance(lengthSums.begin(),which);
        indicies.push_back(distance(lengthSums.begin(),totLen));
        indicies.push_back(index);
        
        std::list<std::list<std::vector<int> > > r;
        for(n = 0; n < (int) baseNrs.size();n++){       // für completeness und contamination 
            
            std::list<std::vector<int> > re;
            accuContigs = 0;
            chromBaseNrs = fromWhichHowMany(minContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n],seed+i);
            justZero = (accumulate(chromBaseNrs.begin(),chromBaseNrs.end(),0) == 0);
            if(justZero){
                Rcout << *((*next(IDs.begin(),indicies[n])).begin()) << " " << i << " " << indicies[n] << " " << partCovered << std::endl;
            }
            if(!justZero){    
                l = 0;
                for(j = 0; j < (int) chromBaseNrs.size();j++){      // für alle chromosomen
                    
                    contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[j],distr,seed+j+n+1);
                    if(contigs.size() <= 0){
                        Rcout << "no Contigs at: " << j << " " << n << " " << i << " " << (*next(IDs.begin(),indicies[n]))[j] << std::endl;
                    }
                    if( contigs.size() > 0){
                        
                        contigSum = accumulate(contigs.begin(),contigs.end(),0);
                        
                        spaces = randomSpaces((int)contigs.size() +1,(*next(lengths.begin(),indicies[n]))[j] -contigSum,seed);
                        co = contigs.begin();
                        sp = spaces.begin();
                        starts.reserve(contigs.size());
                        ends.reserve(contigs.size());
                        
                        swtch = true;
                        at = 0;
                        for(int m = 0;m < (int)contigs.size();m++){     // für alle Contigs
                            at += (*sp);
                            if(swtch){
                                starts.push_back(at+1);
                            }
                            at += (*co);
                            sp++;
                            co++;
                            if((*sp) > 0){
                                ends.push_back(at);
                                swtch = true;
                            }
                            else{
                                swtch = false;
                            }
                        }           // ende alle Contigs
                        
                        at += (*sp);
                        if(ends.size() == 0 || (*sp) == 0){
                            ends.push_back(at);
                        }
                        
                        re.push_back(std::vector<int> {(*next(IDs.begin(),indicies[n]))[j]});
                        re.push_back(starts);
                        re.push_back(ends);
    
                        accuContigs += contigSum;
                        
                        starts.clear();
                        ends.clear();
                    }
                }       // ende alle Chromosomen
                re.push_back(std::vector<int> {accuContigs});
                re.push_back(std::vector<int> {(*totLen)});
                r.push_back(re);
            }
        }       // ende completness und contamination
        if(!justZero){    
            res.push_back(r);
            baseNrs.clear();
            indicies.clear();
        }
    }
    return res;
}