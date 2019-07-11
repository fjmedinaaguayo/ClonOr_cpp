//
//  moveaddedgeMAIS.h
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 01/06/2018.
//  Copyright © 2018 Felipe Medina Aguayo. All rights reserved.
//

#ifndef moveaddedgeMAIS_h
#define moveaddedgeMAIS_h

#include "move.h"

namespace weakarg
{
    
    /**
     @brief This move adds an edge using AISRJ
     */
    class MoveAddEdgeMAIS : public Move
    {
    public:
        MoveAddEdgeMAIS(Param * p,double a);
        Move * clone()
        {
            return new MoveAddEdgeMAIS(*this);
        }
        int move(vector<int> * samplespace=NULL);
        double gammaAIS(int t);
        double logSumExp(vector<double> x);
        vector<int> syst_resamp(vector<double> lw, int N);
        int mult_resamp(vector<double> lw);
        int moveWithinAIS(int t, RecTreeAux* rectreeAux, vector<double> store_ll, double* ll);
        inline int move(){return(move(NULL));}
        ~MoveAddEdgeMAIS();
        
    };
    
} // end namespace weakarg

#endif /* moveaddedgeAIS_h */
