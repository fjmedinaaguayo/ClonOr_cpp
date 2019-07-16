//
//  moveremedgeMAIS.hpp
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 30/04/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#ifndef moveremedgeMAIS_h
#define moveremedgeMAIS_h

#include "move.h"
#include "moveWithinAIS.h"

namespace weakarg
{
    
    /**
     @brief This move adds an edge using AISRJ
     */
    class MoveRemEdgeMAIS : public Move
    {
    public:
        MoveRemEdgeMAIS(Param * p,double a);
        Move * clone()
        {
            return new MoveRemEdgeMAIS(*this);
        }
        int move(vector<int> * samplespace=NULL);
        double logSumExp(vector<double> x);
        vector<int> syst_resamp(vector<double> lw, int N);
        int mult_resamp(vector<double> lw);
        inline int move(){return(move(NULL));}
        ~MoveRemEdgeMAIS();
    };
    
} // end namespace weakarg

#endif
