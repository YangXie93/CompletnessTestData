#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <list>
#include<iterator>
#include <chrono>

using namespace Rcpp;

//[[Rcpp::export]]

List getPfamCountRelation(List refList, std::vector<int> orgVec, std::vector<int> query){
    
    std::vector<int> missing (orgVec.size(),0);
    std::vector<int> exceeding (orgVec.size(),0);
    std::vector<int> total (orgVec.size(),0);
    std::vector<bool> isIn (orgVec.size(),true);
    
    std::vector<int>::iterator missingIt = missing.begin();
    std::vector<int>::iterator exceedingIt = exceeding.begin();
    std::vector<int>::iterator totalIt = total.begin();
    std::vector<bool>::iterator isInIt = isIn.begin();
    
    std::vector<int>::iterator queryIt = query.begin();
    
    List::iterator orgListIt;
    std::vector<int>::iterator pfamCountVecIt;
    std::vector<int>::iterator orgListVecIt;
    
    std::vector<int>::iterator orgIt;
    std::vector<int> countVec;
    std::vector<int> orgV;
    
    
    int dictSize = *std::max_element(orgVec.begin(),orgVec.end());
    std::vector<int> dict (dictSize,-1);
    
    for(int t = 0;t < orgVec.size();t++){
        dict[orgVec[t]] = t;
    }
    
    int i = 0;
    

    for(List::iterator refListIt = refList.begin(); refListIt != refList.end(); refListIt++){
        
        List refListList = as<List>((*refListIt));
        List::iterator orgL_u_countV = refListList.begin();
        
        if(refListList.size() == 0){
            if(*next(queryIt,i) != 0){
                for(int j =0 ;j < exceeding.size();j++){
                    *next(exceedingIt,j) += *next(queryIt,i);
                }
            }
        }
        else{
            countVec = as<std::vector<int> >((*orgL_u_countV));
            orgL_u_countV++;
            List orgL = as<List>(*orgL_u_countV);
            
            pfamCountVecIt = countVec.begin();
            orgListIt = orgL.begin();
            int tmp;
            
            for(int j = 0; j < countVec.size();j++){
                
                orgV = as<std::vector<int> >((*orgListIt));
                
                for(int n = 0; n < orgV.size();n++){
                    
                    tmp = (*pfamCountVecIt) - *next(queryIt,i);
                    
                    int dist = dict[orgV[n]];
                    
                    
                    if(tmp < 0){
                        *next(exceedingIt,dist) += -1* tmp;
                    }
                    else{
                        *next(missingIt,dist) += tmp;
                    }
                    *next(totalIt,dist) += (*pfamCountVecIt);
                }
                
                
                pfamCountVecIt++;
                orgListIt++;
            }
            
        }
       
       
       i++; 
    }
    
    List res = List::create(Named("total") = total,Named("missing") = missing,Named("exceeding") = exceeding);
    
    return res;
}
