//
//  moveremedgeAIS.cpp
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 26/02/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "moveremedgeAIS.h"

using namespace std;
namespace weakarg
{
    
    MoveRemEdgeAIS::MoveRemEdgeAIS(Param * p,double a)
    : Move(p,a)
    {
        description= "Update removing a recombinant edge";
        name= "RemEdge";
        
    }
    
    double MoveRemEdgeAIS::gammaAIS(int t, int T_AIS){
        
        if(t==0)
            return 0;
        
        double result=1-pow(1-(t+0.0)/T_AIS ,5);
        return(result);
    }
    
    int MoveRemEdgeAIS::move(vector<int> * samplespace)
    {
        int T_AIS=param->getT_AIS();
        RecTree * rectree=param->getRecTree();
        int which0=rectree->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
        if (which0<0)
        {
            dlog(1)<<"No valid edges to change!"<<endl;
            return(-1);
        }
        numcalls++;
        
        double lratio=0;
        
        double l0=param->getLL(), l0_original=l0;
        double tfrom0=param->getRecTree()->getRecEdge(which0)->getTimeFrom();
        double tto0  =param->getRecTree()->getRecEdge(which0)->getTimeTo  ();
        int start0=param->getRecTree()->getRecEdge(which0)->getStart();
        int end0  =param->getRecTree()->getRecEdge(which0)->getEnd();
        int efrom0=param->getRecTree()->getRecEdge(which0)->getEdgeFrom();
        int eto0  =param->getRecTree()->getRecEdge(which0)->getEdgeTo  ();
        
        double tfrom_star, tto_star;
        unsigned int start_star, end_star, efrom_star, eto_star;
        int which_star;
        double lu=log(gsl_rng_uniform(rng));
        
        rectree->remRecEdge(which0);
        vector<double> store0(end0-start0);
        for (int i=start0;i<end0;i++)
            store0[i-start0]=param->getLLsite(i);
        
        param->computeLikelihood(start0,end0);
        double l1=param->getLL(), l_star;
        
        
        for(int t=0; t<T_AIS; t++){
            
            if(t==0)
                lratio+=l1-(gammaAIS(t+1,T_AIS)-gammaAIS(t,T_AIS))*l0+log((1.0+rectree->numRecEdge())*2.0/param->getRho()/rectree->getTTotal());
            else
                lratio+=-(gammaAIS(t+1,T_AIS)-gammaAIS(t,T_AIS))*l0;
            
            if(lu<=lratio){
                break;
            }
            
            if(t!=T_AIS-1){
                
                //Draw start and end
                rectree->setBlock(&start_star,&end_star,param->getDelta(),param->getData()->getBlocks());
                //Draw eto and tto
                bool insamplespace=false;
                while(!insamplespace){
                    eto_star=rectree->getPoint(&tto_star);
                    //Draw efrom and tfrom
                    tfrom_star=tto_star+rectree->getNode(eto_star)->getAge();
                    efrom_star=rectree->getEdgeCoal(&tfrom_star);
                    tfrom_star-=rectree->getNode(efrom_star)->getAge();
                    if(samplespace==NULL) insamplespace=true;
                    else if (samplespace->size()==0) insamplespace=true;
                    else
                    {
                        for(unsigned int i=0;i<samplespace->size();i++) if((int)eto_star==samplespace->at(i)||(int)efrom_star==samplespace->at(i))
                        {
                            insamplespace=true;break;
                        }
                    }
                }
                
                vector<double> store_star(end_star-start_star);
                for (unsigned int i=start_star;i<end_star;i++)
                    store_star[i-start_star]=param->getLLsite(i);
                which_star=rectree->addRecEdge(tfrom_star,tto_star,start_star,end_star,efrom_star,eto_star);
                if(which_star<0) return(-1);
                param->computeLikelihood(start_star,end_star);
                l_star=param->getLL();
                
                param->getRecTree()->remRecEdge(which_star);
                for (unsigned int i=start_star;i<end_star;i++)
                    param->setlocLL(i,store_star[i-start_star]);
                
                if (log(gsl_rng_uniform(rng))<=(l_star-l0)*(1-gammaAIS(t+1,T_AIS))){
                    
                    dlog(1)<<"AISRJ t="<<t<<" out of "<<T_AIS<<"..."<<"MCMC in AIS Accepted!"<<endl;
                    l0=l_star;
                    
                }
            }
        }
        
        
        dlog(1)<<"Proposing to remove edge "<<which0<<"...";
        if (lu>lratio)
        {
            dlog(1)<<"Rejected!"<<endl;
            if(param->getRecTree()->addRecEdge(tfrom0,tto0,start0,end0,efrom0,eto0)<0) throw("MoveRemEdge: Can't restore edge!");
            for (int i=start0;i<end0;i++)
                param->setlocLL(i,store0[i-start0]);
            param->setLL(l0_original);
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        param->setLL(l1);
        numaccept++;
        return(1);
    }
    
    MoveRemEdgeAIS::~MoveRemEdgeAIS()
    {}
    //
    
} // end namespace weakarg
