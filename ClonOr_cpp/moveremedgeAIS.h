//
//  moveremedgeAIS.h
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 26/02/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#ifndef MOVEREMEDGEAIS_H
#define MOVEREMEDGEAIS_H

#include "move.h"

namespace weakarg
{
    
    /**
     @brief This move removes an edge
     */
    class MoveRemEdgeAIS : public Move
    {
    public:
        MoveRemEdgeAIS(Param * p,double a);
        Move * clone()
        {
            return new MoveRemEdgeAIS(*this);
        }
        int move(vector<int> * samplespace=NULL);
        double gammaAIS(int t, int up);
        inline int move(){return(move(NULL));}
        ~MoveRemEdgeAIS();
        
    };
    
} // end namespace weakarg
#endif
