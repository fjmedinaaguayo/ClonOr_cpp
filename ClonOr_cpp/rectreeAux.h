//
//  rectreeAux.h
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 11/06/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#ifndef rectreeAux_h
#define rectreeAux_h

#include <stdio.h>
#include "rectree.h"
#include <cmath>
#define LOG025 -1.38629461

namespace weakarg
{
    
class RecTreeAux{
public:
    Data* data;
    RecTree* rectree0;
    double tfrom, tto;
    unsigned int start, end, efrom, eto;
    int which;
    vector<double> store_ll;
    std::vector<int> tabNode;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabRec;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabFather;///<Used only in makeLocalTree (made global for speed)
    std::vector<double> age;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabEdge;///<Used only in makeLocalTree (made global for speed)
    std::vector<std::vector<int> > tabSons;
    std::vector<std::vector<double> > tabSonsDist;
    
    
    RecTreeAux(Data* data, RecTree* rectree0);
    RecTreeAux(RecTreeAux* Copy);
    ~RecTreeAux();
    void setAll(unsigned int start, unsigned int end, double tfrom, double tto, unsigned int efrom, unsigned int eto);
    double computePartialLL(double theta);
    void computeLocalTree(double theta, unsigned int site, std::vector<RecEdge*> edge);
};

}

#endif /* rectreeAux_h */
