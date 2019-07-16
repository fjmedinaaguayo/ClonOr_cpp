//
//  moveWithinAIS.hpp
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 11/07/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#ifndef moveWithinAIS_hpp
#define moveWithinAIS_hpp

#include "param.h"
#include "rectreeAux.h"


namespace weakarg
{
    class MoveWithinAIS
    {
        public:
            MoveWithinAIS(int T, double pow);
            double gammaAIS(int t);
            int moveFwd(int t, Param* param, RecTreeAux* rectreeAux, double* ll);
            int moveStartEnd(int t, Param* param, RecTreeAux* rectreeAux, double* ll);
            int moveTimes(int t, Param* param, RecTreeAux* rectreeAux, double* ll);
            int moveFull(int t, Param* param, RecTreeAux* rectreeAux, double* ll);
            int movePrior(int t, Param* param, RecTreeAux* rectreeAux, double* ll);
        
        private:
            int AIS_T;
            double AIS_power;
        
    };
}

#endif
