//
//  moveaddedgeAIS.cpp
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 01/06/2018.
//  Copyright © 2018 Felipe Medina Aguayo. All rights reserved.
//

#include "moveaddedgeAIS.h"
//
using namespace std;
namespace weakarg
{
    
    MoveAddEdgeAIS::MoveAddEdgeAIS(Param *p,double a)
    : Move(p,a)
    {
        description= "Update adding a recombinant edge using AISRJ";
        name= "AddEdge";
    }
    
    double MoveAddEdgeAIS::gammaAIS(int t){
        
        if(t==0)
            return 0;
        
        double result=pow((t+0.0)/param->getT_AIS() , param->getgamma_AIS());
        return(result);
    }
    
    int MoveAddEdgeAIS::move(vector<int> * samplespace)
    {
        int T_AIS=param->getT_AIS();
        RecTree * rectree=param->getRecTree();
        
        double tfrom1,tto1;
        double tfrom2,tto2;
        unsigned int start1,end1;
        unsigned int start2,end2;
        unsigned int efrom1,eto1;
        unsigned int efrom2,eto2;
        int which1;
        int which2;
        double lratio=0, l0, l1, l2;
        
        //Draw start and end
        rectree->setBlock(&start1,&end1,param->getDelta(),param->getData()->getBlocks());
        //Draw eto and tto
        bool insamplespace=false;
        while(!insamplespace){
            eto1=rectree->getPoint(&tto1);
            //Draw efrom and tfrom
            tfrom1=tto1+rectree->getNode(eto1)->getAge();
            efrom1=rectree->getEdgeCoal(&tfrom1);
            tfrom1-=rectree->getNode(efrom1)->getAge();
            if(samplespace==NULL) insamplespace=true;
            else if (samplespace->size()==0) insamplespace=true;
            else
            {
                for(unsigned int i=0;i<samplespace->size();i++) if((int)eto1==samplespace->at(i)||(int)efrom1==samplespace->at(i))
                {
                    insamplespace=true;break;
                }
            }
        }
        
        l0=param->getLL();
        vector<double> store1(end1-start1);
        for (unsigned int i=start1;i<end1;i++)
            store1[i-start1]=param->getLLsite(i);
        which1=rectree->addRecEdge(tfrom1,tto1,start1,end1,efrom1,eto1);
        if(which1<0) return(-1);
        param->computeLikelihood(start1,end1);
        double lu=log(gsl_rng_uniform(rng));
        
        for(int t=0; t<T_AIS; t++){
            
            l1=param->getLL();
            
            if(t==0)
                lratio+=(gammaAIS(t+1)-gammaAIS(t))*l1-l0+log(param->getRho()*rectree->getTTotal()/2.0/rectree->numRecEdge());
            else
                lratio+=(gammaAIS(t+1)-gammaAIS(t))*l1;
            
            if(lratio<lu){
                break;
            }
            
            if(t!=T_AIS-1){
                
                param->getRecTree()->remRecEdge(which1);
                for (unsigned int i=start1;i<end1;i++)
                    param->setlocLL(i,store1[i-start1]);
                param->setLL(l0);
                
                //Draw start and end
                rectree->setBlock(&start2,&end2,param->getDelta(),param->getData()->getBlocks());
                //Draw eto and tto
                bool insamplespace=false;
                while(!insamplespace){
                    eto2=rectree->getPoint(&tto2);
                    //Draw efrom and tfrom
                    tfrom2=tto2+rectree->getNode(eto2)->getAge();
                    efrom2=rectree->getEdgeCoal(&tfrom2);
                    tfrom2-=rectree->getNode(efrom2)->getAge();
                    if(samplespace==NULL) insamplespace=true;
                    else if (samplespace->size()==0) insamplespace=true;
                    else
                    {
                        for(unsigned int i=0;i<samplespace->size();i++) if((int)eto2==samplespace->at(i)||(int)efrom2==samplespace->at(i))
                        {
                            insamplespace=true;break;
                        }
                    }
                }
                
                vector<double> store2(end2-start2);
                for (unsigned int i=start2;i<end2;i++)
                    store2[i-start2]=param->getLLsite(i);
                which2=rectree->addRecEdge(tfrom2,tto2,start2,end2,efrom2,eto2);
                if(which2<0) return(-1);
                param->computeLikelihood(start2,end2);
                l2=param->getLL();
                
                if (log(gsl_rng_uniform(rng))>(l2-l1)*(gammaAIS(t+1))){
                    
                    param->getRecTree()->remRecEdge(which2);
                    for (unsigned int i=start2;i<end2;i++)
                        param->setlocLL(i,store2[i-start2]);
                    param->setLL(l0);
                    
                    which1=rectree->addRecEdge(tfrom1,tto1,start1,end1,efrom1,eto1);
                    if(which1<0) return(-1);
                    param->computeLikelihood(start1,end1);
                }
                else{
                    
                    dlog(1)<<"AISRJ t="<<t<<" out of "<<T_AIS<<"..."<<"MCMC in AIS Accepted!"<<endl;
                    tfrom1=tfrom2;
                    tto1=tto2;
                    start1=start2;
                    end1=end2;
                    efrom1=efrom2;
                    eto1=eto2;
                    l1=l2;
                    which1=which2;
                    
                    store1=store2;
                }
            }
        }
        
        numcalls++;
        dlog(1)<<"Proposing to add edge via AISRJ "<<efrom1<<":"<<tfrom1<<"->"<<eto1<<":"<<tto1<<"...";
        
        if (lu>lratio)
        {
            dlog(1)<<"Rejected!"<<endl;
            param->getRecTree()->remRecEdge(which1);
            for (unsigned int i=start1;i<end1;i++)
                param->setlocLL(i,store1[i-start1]);
            param->setLL(l0);
            return(0);
            //param->computeLikelihood(start,end);
        }
        else dlog(1)<<"Accepted!"<<endl;
        numaccept++;
        return(1);
    }
    
    MoveAddEdgeAIS::~MoveAddEdgeAIS()
    {}
    //
    
} // end namespace weakarg
