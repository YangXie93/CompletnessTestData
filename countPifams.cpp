#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <Rcpp.h>
#include <time.h>

using namespace std;


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

Rcpp::List countPifams(list<list<vector<int> > >& pifams,list<list<vector<int> > >& ORFs,list<list<list<vector<int> > > >& contigs,vector<int>& names){
    Rcpp::Rcout << "countPifams begin" << endl;
    Rcpp::Rcout << "-----------------" << endl;

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

    vector<int> tmp1;

    vector<int>::iterator trans;
    vector<int>::iterator t;
    vector<int>::iterator starts;
    vector<int>::iterator ends;
    vector<int>::iterator width;
    vector<int>::iterator values;
    vector<int>::iterator orfValues;
    vector<int>::iterator orfWidths;


    list<vector<int> >::iterator orfs;
    list<vector<int> >::iterator pis;

    list<list<list<vector<int> > > >::iterator con;
    list<list<vector<int> > >::iterator cons;
    list<vector<int> >::iterator conts;

    vector<int>::iterator bN;
    vector<bool>::iterator nio;
    vector<int>::iterator tms;


    for(con = contigs.begin();con != contigs.end();con++){
        // alle Werte
        Rcpp::List tmp;

        for(cons = (*con).begin();cons != (*con).end();cons++){
            // alle Genome
            for(conts = (*cons).begin();conts != (*cons).end();conts = next(conts,3)){
                // alle Chromosomen

                name.push_back(*((*conts).begin()));
                trans = find(names.begin(),names.end(),*((*conts).begin()));
                if(trans != names.end()){
                    transfer = distance(names.begin(),trans);
                    orfs = (*next(ORFs.begin(),transfer)).begin();
                    for(pis = (*next(pifams.begin(),transfer)).begin();pis != (*next(pifams.begin(),transfer)).end();pis = next(pis,2)){
                        // für alle 6 GENOME und ORF Werte

                        starts = (*next(conts)).begin();
                        ends = (*next(conts,2)).begin();
                        width = (*pis).begin();
                        values = (*next(pis)).begin();
                        orfValues = (*next(orfs)).begin();

                        for(orfWidths = (*orfs).begin(); distance(orfWidths,(*orfs).end()) > 0; orfWidths++){
                            //für alle Length Werte in ORF
                            if(starts != (*next(conts)).end() && (*orfValues) >= 0 && ((*starts) < o+(*orfWidths)) ){
                                fill(notInOrf.begin(),notInOrf.end(),true);

                                while(width != (*pis).end() && o+(*orfWidths) > n){


                                    if(n >= o && (*values) >= 0){
                                        while(starts != (*next(conts)).end() && n+(*width) > (*starts)){

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
                                            if((*ends) >= (n+(*width))){
                                                break;
                                            }
                                            else{
                                                x += (*ends) -(*starts) +1;
                                                starts++;
                                                ends++;
                                            }
                                        }

                                    }
                                    n += (*width);
                                    values++;
                                    width++;
                                }

                            }
                            o += (*orfWidths);
                            orfValues++;
                        }
                        //end alle Length Werte aus ORF
                        if(pis == (*next(pifams.begin(),transfer)).begin()){
                            total += o;
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
                    if(cons == (*con).begin()){
                        partial += zwsch;
                    }
                    else{
                        contTotal += zwsch;
                    }
                    zwsch = 0;
                }
            }
            //ende alle Chromosomen
            if(cons == (*con).begin()){
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
        }
        //ende alle Genome
        res.push_back(tmp);
    }
    Rcpp::Rcout << "countPifams end" << endl;
    return res;
}
