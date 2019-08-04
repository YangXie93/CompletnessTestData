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

Rcpp::List countPifams(list<list<vector<int> > >& pifams,list<list<list<vector<int> > > >& contigs){

    Rcpp::List res;
    Rcpp::IntegerVector id;
    Rcpp::IntegerVector times;
    Rcpp::IntegerVector baseNum;
    int n = 0;
    for(list<list<list<vector<int> > > >::iterator con = contigs.begin();con != contigs.end();con++){
        Rcpp::List tmp;
        for(list<list<vector<int> > >::iterator cons = (*con).begin();cons != (*con).end();cons++){
            Rcpp::List tmp2;
            for(list<vector<int> >::iterator conts = (*cons).begin();conts != (*cons).end();conts = next(conts,2)){
                list<list<vector<int> > >::iterator pi = pifams.begin();
                int which = (*((*conts).begin()));
                for(list<vector<int> >::iterator pis = (*next(pi,which)).begin();pis != (*next(pi,which)).end();pis = next(pis,2)){
                    vector<int>::iterator starts = (*next(conts)).begin();
                    vector<int>::iterator ends = (*next(conts,2)).begin();
                    vector<int>::iterator values = (*next(pis)).begin();
                    
                    for(vector<int>::iterator width = (*pis).begin();width != (*pis).end();width++){
                        if((*values) >= 0){
                            vector<int> tmp = intervallOverlap(n,n+(*width),(*starts),(*ends));
                            if((int) tmp.size() > 1){
                                Rcpp::IntegerVector::iterator t = find(id.begin(),id.end(),(*values));
                                if(t == id.end()){
                                    baseNum.push_back(tmp[1] - tmp[0]);
                                    times.push_back(1);
                                    id.push_back((*values));
                                }
                                else{
                                    times[((int) (t-id.begin()))]++;
                                    baseNum[((int) (t-id.begin()))] += tmp[1] - tmp[0];
                                }
                            }
                        }
                        n += (*width);
                        if(n >= (*ends)){
                            starts++;
                            ends++;
                        }
                        values++;
                    }
                    n = 0;
                }
                if(cons != (*con).begin()){
                    tmp2.push_back(Rcpp::List::create(Rcpp::Named("compPifamNames") = id,Rcpp::Named("compPifamCount") = times,Rcpp::Named("compBaseCount") = baseNum));
                }
                else{
                    tmp2.push_back(Rcpp::List::create(Rcpp::Named("contPifamNames") = id,Rcpp::Named("contPifamCount") = times,Rcpp::Named("contBaseCount") = baseNum));
                }
            }
            tmp.push_back(tmp2);
        }
        con++;
        res.push_back(tmp);
    }
    return res;
}




Rcpp::List contPifams(list<list<vector<int> > >& pifams,list<list<list<vector<int> > > >& contigs){
    
    Rcpp::List res;
    Rcpp::IntegerVector id;
    Rcpp::IntegerVector times;
    Rcpp::IntegerVector baseNum;
    for(list<list<list<vector<int> > > >::iterator con = contigs.begin();con != contigs.end();con++){
        for(list<list<vector<int> > >::iterator cons = (*con).begin();cons != (*con).end();cons++){
            for(list<vector<int> >::iterator cont = (*cons).begin();cont != (*cons).end();cont++){
                
            }
        }
    }
    
    return res;
}


