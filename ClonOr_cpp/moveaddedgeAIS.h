//
//  moveaddedgeAIS.h
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 01/06/2018.
//  Copyright © 2018 Felipe Medina Aguayo. All rights reserved.
//

#ifndef moveaddedgeAIS_h
#define moveaddedgeAIS_h

#include "move.h"

namespace weakarg
{
    
    /**
     @brief This move adds an edge using AISRJ
     */
    class MoveAddEdgeAIS : public Move
    {
    public:
        MoveAddEdgeAIS(Param * p,double a);
        Move * clone()
        {
            return new MoveAddEdgeAIS(*this);
        }
        int move(vector<int> * samplespace=NULL);
        double gammaAIS(int t);
        inline int move(){return(move(NULL));}
        ~MoveAddEdgeAIS();
        
    };
    
} // end namespace weakarg

#endif /* moveaddedgeAIS_h */
