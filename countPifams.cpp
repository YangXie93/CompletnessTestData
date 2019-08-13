#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <Rcpp.h>

using namespace std;
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
    
    for(list<list<list<vector<int> > > >::iterator con = contigs.begin();con != contigs.end();con++){
       
        Rcpp::List tmp;
        int total = 0; 
        int partial = 0;
        int contTotal = 0;
        
        for(list<list<vector<int> > >::iterator cons = (*con).begin();cons != (*con).end();cons++){
            vector<int> name;
            
            for(list<vector<int> >::iterator conts = (*cons).begin();conts != (*cons).end();conts = next(conts,3)){
                
                int o = 1;
                int n = 1;
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
                        t += (*e)-(*s);
                        e++;
                    }
                    if(conts == (*cons).begin()){
                        partial = t;
                    }
                    contTotal += partial;
                    Rcpp::Rcout << "Contig Length: " << contTotal << endl;
                    for(list<vector<int> >::iterator pis = (*next(pifams.begin(),transfer)).begin();pis != (*next(pifams.begin(),transfer)).end();pis = next(pis,2)){
                        // für alle 6 GENOME und ORF Werte
                        
                        list<vector<int> >::iterator orfs = (*next(ORFs.begin(),transfer)).begin();
                        vector<int>::iterator starts = (*next(conts)).begin();
                        vector<int>::iterator ends = (*next(conts,2)).begin();
                        vector<int>::iterator width = (*pis).begin();
                        vector<int>::iterator values = (*next(pis)).begin();
                        vector<int>::iterator orfValues = (*next(orfs)).begin();
                        int test = 0;
                        for(vector<int>::iterator orfWidths = (*orfs).begin(); distance(orfWidths,(*orfs).end()) > 0; orfWidths++){
                            //für alle Length Werte in ORF
                            if(starts != (*next(conts)).end()){
                                while(width != (*pis).end() && orfWidths != (*orfs).end() && o+(*orfWidths) > n){
                                    
                                    
                                    if((*orfValues) >= 0){
                                        for(int i = 0;i < (int) notInOrf.size();i++){
                                            notInOrf[i] = true;
                                        }
                                        while(starts != (*next(conts)).end() && width != (*pis).end() && n+(*width) > (*starts)){
                                            
                                            if((*values) >= 0){
                                                vector<int> tmp = intervallOverlap(n,n+(*width),(*starts),(*ends));
                                                if((int) tmp.size() > 1){
                                                    vector<int>::iterator t = find(id.begin(),id.end(),(*values));
                                                    if(t == id.end()){
                                                        test++;
                                                        baseNum.push_back(tmp[1] - tmp[0]);
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
                                            x += (*ends) -(*starts);
                                            starts++;
                                            ends++;
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
                        Rcpp::Rcout << "X: " << x << " n: " << n << " o: " << o << endl;
                        Rcpp::Rcout << "res size: " << id.size() << " should be: "<< test << endl;
                        o = 1;
                        n = 1;
                        x = 0;
                        orfs = next(orfs,2);
                    }
                    //ende alle 6 GENOME und ORF Werte
                    
                }
            }
            
            Rcpp::Rcout << "---------------------------------------------" << endl;
            if(cons == (*con).begin()){
                vector<double> completeness = {(double)partial/(double) total};
                tmp = (Rcpp::List::create(Rcpp::Named("compChromID") = name,Rcpp::Named("compPifamNames") = id,Rcpp::Named("compPifamCount") = times,Rcpp::Named("compBaseCount") = baseNum,Rcpp::Named("completness") = completeness));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
            else{
                vector<double> contamination = {((double)partial/(double)contTotal)};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("contChromID") = name,Rcpp::Named("contPifamNames") = id,Rcpp::Named("contPifamCount") = times,Rcpp::Named("contBaseCount") = baseNum,Rcpp::Named("contamination") = contamination));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
                notInOrf.clear();
            }
        }
        total = 0;
        contTotal = 0;
        partial = 0;
        res.push_back(tmp);
    }
    return res;
}

