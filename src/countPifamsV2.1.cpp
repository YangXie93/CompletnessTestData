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
        return std::vector<int> (1,0);
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




std::vector<int> fromWhichHowMany(int minContigLength,int totalLength,std::vector<int> lengths,int needed,int seed){
    std::default_random_engine generator;
    generator.seed(seed);
    std::vector<int> res;
    res.reserve(lengths.size()/2);
    int share;
    int tmp;
    std::vector<int> tmp1;
    std::vector<int> accession;
    tmp1.reserve(lengths.size()/2);
    accession.reserve(lengths.size()/2);
    int j = 0;
    for(int i = 1;i < (int) lengths.size();i+= 2){
        res.push_back(0);
        if(lengths[i] >= minContigLength){
            tmp1.push_back(lengths[i]);
            accession.push_back(j);
        }
        else{
            totalLength -= lengths[i];
        }
        j++;
    }
    if(tmp1.size() > 0){
        std::vector<unsigned int> tmp2;
        for(int i =0; i < tmp1.size();i++){
            tmp2.push_back(tmp1[i] *(generator()% 100));
        }
        unsigned int ges = accumulate(tmp2.begin(),tmp2.end(),0);
        for(int i = 0;i < tmp2.size();i++ ){
            res[accession[i]] = needed *(tmp2[i]/(double) ges);
        }
    }
    return res;
}

//' @export
//[[Rcpp::export]]
std::vector<std::vector<int> > mkContigs(std::list<std::vector<int> >& lengths,std::vector<int>& lengthSums,std::vector<int>& total,std::vector<int>& partial,int minContigLength,int meanContigLength,int number,std::vector<double>& comp,std::vector<double>& cont,std::vector<std::vector<int> >& names,std::vector<int>& access,int seed = 0,std::string distr = "normal"){

    std::default_random_engine generator;
    generator.seed(seed);

    std::vector<std::vector<int> > res;
    res.reserve(number);

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
    std::vector<int>::iterator nms;

    int at = 0;
    double partCovered;
    std::vector<int>::iterator which;
    std::vector<int>::iterator max;
    bool swtch;
    int index;
    int l;
    int n;
    int j;
    int count;
    int accuContigs;

    if(cont[1] > comp[0]){
        cont[1] = comp[0];
    }
    if(cont[0] == cont[1]){
        cont[0] = cont[1] -0.01;
    }
    if(comp[0] == comp[1]){
        comp[0] = comp[1] -0.01;
    }
    for(int i = 0; i < number; i++){

        count = 0;
        
        partCovered = ((comp[0] *100) + (generator() % (int)((comp[1]-comp[0]) *100)))/100.0;
        which = next(lengthSums.begin(),(generator() % lengthSums.size()));
        if(i == 6501){
            Rcout << partCovered << std::endl;
        }
        
        while(((*which)* partCovered) < minContigLength && count < (int)lengthSums.size()){
            which++;
            count++;
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
        double contPart = ((cont[0] *100) + (generator() % ((int)(cont[1]*100) - (int)(cont[0] *100))))/100.0;
        if((*totLen) *contPart <= minContigLength){
            if(contPart != 0){
                baseNrs.push_back(minContigLength);
            }
        }
        else{
            baseNrs.push_back((*totLen) *contPart);
        }

        
        count = 0;
        which = next(lengthSums.begin(),(generator() % (int)lengthSums.size()));
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
            }
        }
        if(((*which) < *prev(baseNrs.end()) && count == (int)lengthSums.size()) || which == lengthSums.end()){
            which = max;
            baseNrs.pop_back();
            baseNrs.push_back((*which) *0.9);
        }
        
        
        res.push_back(std::vector<int> (0));
        index = distance(lengthSums.begin(),which);
        indicies.push_back(distance(lengthSums.begin(),totLen));
        indicies.push_back(index);
        
        total.push_back((*totLen));
        
        for(n = 0; n < (int) baseNrs.size();n++){
            
            accuContigs = 0;
            
            res.push_back(std::vector<int> (0));
            chromBaseNrs = fromWhichHowMany(minContigLength,(*next(lengthSums.begin(),indicies[n])),(*next(lengths.begin(),indicies[n])),baseNrs[n],seed);
            if(i == 6501){
                Rcout << "fwhm:" << std::endl;
                for(int p = 0;p < chromBaseNrs.size();p++){
                    Rcout << chromBaseNrs[p] << " : " << (*next(lengths.begin(),indicies[n]))[p*2+1] << std::endl;
                }
                Rcout << std::endl << std::endl;
            }
            nms = (*next(names.begin(),indicies[n])).begin();
            l = 0;
            for(j = 0; j < (int) chromBaseNrs.size();j++){
                if(chromBaseNrs[j] == (*next(lengths.begin(),indicies[n]))[j]){
                    Rcout << i << ": " << chromBaseNrs[j] << " index: " << indicies[n] << "\n";
                }
                contigs = randomContigs(minContigLength,meanContigLength,chromBaseNrs[j],distr,seed+j+n+1);
                if( contigs.size() > 0 && contigs[0] != 0){
                    int contigSum = accumulate(contigs.begin(),contigs.end(),0);
                    accuContigs += contigSum;
                    spaces = randomSpaces((int)contigs.size() +1,(*next(lengths.begin(),indicies[n]))[j *2+1] -contigSum,seed);
                    co = contigs.begin();
                    sp = spaces.begin();
                    starts.reserve(contigs.size());
                    ends.reserve(contigs.size());
                    swtch = true;
                    at = 0;
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
                                swtch = true;
                            }
                            else{
                              swtch = false;
                            }
                    }
                    at += (*sp);
                    if(ends.size() == 0 || (*sp) == 0){
                        ends.push_back(at);
                    }

                    res.push_back(std::vector<int> {(*next(lengths.begin(),indicies[n]))[j *2]});
                    res.push_back(starts);
                    res.push_back(ends);
                    
                    access.push_back(*next(nms,l));
                    starts.clear();
                    ends.clear();
                }
            }
            partial.push_back(accuContigs);
        }
        baseNrs.clear();
        indicies.clear();
    }
    return res;
}




int intervallOverlap(int start1, int end1,int start2,int end2){

    if(start2 > end1 || end2 < start1 ){
        return 0;
    }
    if(start1 > start2){
        if(end2 < end1){
            return end2 -start1 +1;
        }
        else{
            return end1 - start1+1;
        }
    }
    else{
        if(end2 < end1){
            return end2 - start2+1;
        }
        else{
            return end1 - start2+1;
        }
    }
}

//' @export
//[[Rcpp::export]]

List countPifams(std::list<std::list<std::vector<int> > > &pifams,std::list<std::list<std::vector<int> > > &ORFs,std::vector<std::vector<int> > &contigs,std::vector<int> &names,std::vector<int> &total,std::vector<int>& partial){


    List res;
    std::vector<int> id;
    std::vector<int> times;
    std::vector<int> baseNum;
    std::vector<bool> notInOrf;
    std::vector<int> name;
    std::vector<double> completeness;
    std::vector<double> contamination;
    std::vector<int>::iterator acc = names.begin();
    std::vector<int>::iterator tot = total.begin();
    std::vector<int>::iterator part = partial.begin();

    long o = 0;
    long n = 0;
    int transfer;
    int GenNr;

    int tmp1;

    std::vector<int>::iterator trans;
    std::vector<int>::iterator t;
    std::vector<int>::iterator starts;
    std::vector<int>::iterator ends;
    std::vector<int>::iterator width;
    std::vector<int>::iterator values;
    std::vector<int>::iterator orfValues;
    std::vector<int>::iterator orfWidths;

    std::list<std::vector<int> >::iterator orfs;
    std::list<std::vector<int> >::iterator pis;

    std::vector<std::vector<int> >::iterator con = contigs.begin();

    std::vector<int>::iterator bN;
    std::vector<bool>::iterator nio;
    std::vector<int>::iterator tms;
    bool swtch = false;

    while(distance(con,contigs.end()) > 0){         // alle Werte
        
        List tmp;
        swtch = false;
        GenNr = 1;
        
        while((int)(*con).size() != 0 && distance(con,contigs.end()) > 0){          // alle Genome
            
            swtch = true;
            
            while((int) (*con).size() != 0 && distance(con,contigs.end()) > 0){         // alle Chromosomen
                
                transfer = (*acc);
                name.push_back(*((*con).begin()));
                
                if(transfer != (int)pifams.size()){
                    
                    orfs = (*next(ORFs.begin(),transfer)).begin();
                    
                    for(pis = (*next(pifams.begin(),transfer)).begin();distance(pis,(*next(pifams.begin(),transfer)).end()) > 0;pis = next(pis,2)){     // für alle 6 GENOME und ORF Werte
                        
                        starts = (*next(con)).begin();
                        ends = (*next(con,2)).begin();
                        width = (*pis).begin();
                        values = (*next(pis)).begin();
                        orfWidths = (*orfs).begin();
                        orfValues = (*next(orfs)).begin();
                        o = 0;
                        n = 0;

                        for(starts = (*next(con)).begin();starts != (*next(con)).end();starts++){       //für alle Contigs
                            
                            
                            while((*starts) > (o +(*orfWidths)) && orfWidths != (*orfs).end()){         // weitersetzen des ORF intervalls bis es den Contig schneidet
                                o += (*orfWidths);
                                orfWidths++;
                                orfValues++;
                            }
                            
                            while(o <  (*ends) &&  orfWidths != (*orfs).end()){         //weitersetzen des ORF intervalls solange es den Contig schneidet
                                
                                if((*orfValues) > 0){
                                    
                                    fill(notInOrf.begin(),notInOrf.end(),true);
                                   
                                    while(o > n+(*width) && width != (*pis).end()){         //weitersetzen des pfam intervalls bis es den Contig schneidet
                                        n += (*width);
                                        width++;
                                        values++;
                                    }
                                    
                                    while(n < (*ends) && width != (*pis).end()){         //weitersetzen des pfam intervalls solange es den Contig schneidet
                                        
                                        if((*values) > 0){                                                  //wenn nicht -1
                                            tmp1 = intervallOverlap(n+1,n+(*width),(*starts),(*ends));      //schnitt berechnen
                                            

                                            //--------------------------------------------------- Werte eintragen ------------------------
                                            if(tmp1 > 0){
                                                t = find(id.begin(),id.end(),(*values));
                                                
                                                if(t == id.end()){
                                                    baseNum.push_back(tmp1);
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
                                                    (*bN) += tmp1;
                                                }
                                            }
                                            //---------------------------------------------------------------------------------------------
                                        }
                                        
                                        if((n + (*width)) < (*ends)){
                                            n += (*width);
                                            width++;
                                            values++;
                                            while(n+(*width) >= o && orfWidths != (*orfs).end()){
                                                o += (*orfWidths);
                                                orfWidths++;
                                                orfValues++;
                                            }
                                        }
                                        else{
                                            break;
                                        }
                                    }
                                }
                                
                                if((o + (*orfWidths)) < (*ends) && orfWidths != (*orfs).end()){
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
                        
                        o = 0;
                        n = 0;
                        orfs = next(orfs,2);
                    }
                    //ende alle 6 GENOME und ORF Werte

                }
                acc++;
                con = next(con,3);
                
            }
            //ende alle Chromosomen
            if(GenNr == 1){
                completeness = {(*part)/(double) (*tot)};
                tmp.push_back(List::create(Named("compChromID") = name,Named("compPifamNames") = id,Named("compPifamCount") = times,Named("compBaseCount") = baseNum,Named("completness") = completeness));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
            else{
                contamination = {(*part)/(double)(*tot)};
                tmp.push_back(List::create(Named("contChromID") = name,Named("contPifamNames") = id,Named("contPifamCount") = times,Named("contBaseCount") = baseNum,Named("contamination") = contamination));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
            con++;
            GenNr++;
            part++;
        }
        //ende alle Genome
        if(swtch){
            res.push_back(tmp);
            tot++;
        }
        con++;
    }

    return res;
}


//[[Rcpp::export]]

List compTestData(std::list<std::list<std::vector<int> > > &pifams,std::list<std::list<std::vector<int> > > &ORFs,std::list<std::vector<int> > &lengths,std::vector<int> &lengthSums,int minContigLength,int meanContigLength,int number,std::vector<double> &comp,std::vector<double> &con,std::vector<std::vector<int> > &names,int seed = 0,std::string distr = "normal"){
    
    
    std::vector<int> access;        // indices zu den Daten in ORFs und pifams zu jedem erzeugten Beispiel in mkContigs
    
    std::vector<int> total;         // die lengthSums Werte zu jedem erzeugten beispiel in mkContigs
    
    std::vector<int> partial;       // die ContigSummen zu jedem erzeugtem Beispiel in mkContigs
    
    std::vector<std::vector<int> > conts = mkContigs(lengths,lengthSums,total,partial,minContigLength,meanContigLength,number,comp,con,names,access,seed,distr);
    
    return countPifams(pifams,ORFs,conts,access,total,partial);
}


//' @export
//[[Rcpp::export]]
std::vector<int> chooseGenomesCpp(std::vector<int> lengthSums,int minContigLength,int seed = 1,int times = 100){

    std::vector<int> res;
    std::default_random_engine generator;
    generator.seed(seed);
    std::vector<double> comp = {0.6,1};
    for(int i = 0;i < times;i++){

        int partCovered = ((comp[0] *100) + (generator() % (int)((comp[1]-comp[0]) *100)))/100.0;
        std::vector<int>::iterator which = next(lengthSums.begin(),(generator() % lengthSums.size()));
        int count = 0;

        while(((*which)* partCovered) < minContigLength && count < (int)lengthSums.size()){
            which++;
            count++;
            if(which == lengthSums.end()){
                which = lengthSums.begin();
            }
        }

        res.push_back((distance(lengthSums.begin(),which) +1));
    }

    return res;
}

//' @export
//[[Rcpp::export]]

int testDefaultRandomEngine(int seed){
    std::default_random_engine generator;
    generator.seed(seed);
    return generator();
}
