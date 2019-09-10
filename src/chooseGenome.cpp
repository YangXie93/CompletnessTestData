#include <Rcpp.h>
#include <vector>
#include <random>


//[[Rcpp::export]]
int return10(){
    return 10;
}

//[[Rcpp::export]]

std::vector<int> chooseGenomes(std::vector<int> lengthSums,int minContigLength,int seed = 1,int times = 100){
    
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