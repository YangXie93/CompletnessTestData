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
                int transfer = distance(names.begin(),trans);
                list<vector<int> >::iterator pis = (*next(pifams.begin(),transfer)).begin();
                list<vector<int> >::iterator orfs = (*next(ORFs.begin(),transfer)).begin();
                vector<int>::iterator starts = (*next(conts)).begin();
                vector<int>::iterator ends = (*next(conts,2)).begin();
                vector<int>::iterator width = (*pis).begin();
                vector<int>::iterator values = (*next(pis)).begin();
                vector<int>::iterator orfValues = (*next(orfs)).begin();
                

                for(vector<int>::iterator orfWidths = (*orfs).begin(); orfWidths != (*orfs).end();orfWidths++){
                    if(starts != (*next(conts)).end()){
                        while(o+(*orfWidths) > n && width != (*pis).end()){
                            if((*orfValues) >= 0){
                                vector<bool> notInOrf(true,id.size());
                                while(n+(*width) > (*starts) && starts != (*next(conts)).end()){
                                    if((*values) >= 0){
                                        vector<int> tmp = intervallOverlap(n,n+(*width),(*starts),(*ends));
                                        if((int) tmp.size() > 1){
                                            vector<int>::iterator t = find(id.begin(),id.end(),(*values));
                                            if(t == id.end()){
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
                                    if(n+(*width) >= (*ends)){
                                        x += (*ends) -(*starts);
                                        starts++;
                                        ends++;
                                    }
                                    else{
                                        n+= (*width);
                                        width++;
                                        values++;
                                    }
                                }
                            }
                            if(o+(*orfWidths) >= n+(*width) ){
                                n += (*width);
                                values++;
                                width++;
                            }
                            else{
                                o+= (*orfWidths);
                                orfWidths++;
                                orfValues++;
                            }
                        }
                    o += (*orfWidths);
                    orfValues++;
                    }
                }
                if(cons == (*con).begin()){
                    total += n;
                    partial += x;
                }
                contTotal += x;
                o = 1;
                n = 1;
                x = 0;
            }
            if(cons == (*con).begin()){
                vector<double> completeness = {(double)partial/(double) total};
                tmp = (Rcpp::List::create(Rcpp::Named("compChromID") = name,Rcpp::Named("compPifamNames") = id,Rcpp::Named("compPifamCount") = times,Rcpp::Named("compBaseCount") = baseNum,Rcpp::Named("completness") = completeness));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
            }
            else{
                vector<double> contamination = {(1- (double)partial/(double)contTotal)};
                tmp.push_back(Rcpp::List::create(Rcpp::Named("contChromID") = name,Rcpp::Named("contPifamNames") = id,Rcpp::Named("contPifamCount") = times,Rcpp::Named("contBaseCount") = baseNum,Rcpp::Named("contamination") = contamination));
                name.clear();
                id.clear();
                baseNum.clear();
                times.clear();
            }
        }
        total = 0;
        contTotal = 0;
        partial = 0;
        res.push_back(tmp);
    }
    return res;
}

