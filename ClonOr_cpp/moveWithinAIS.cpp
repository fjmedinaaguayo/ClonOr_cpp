//
//  moveWithinAIS.cpp
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 11/07/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "moveWithinAIS.h"
#include <cmath>
#include "gsl/gsl_sf.h"

using namespace std;
namespace weakarg
{
    MoveWithinAIS::MoveWithinAIS(int T, double pow): AIS_T(T), AIS_power(pow){
        
    }
    
    double MoveWithinAIS::gammaAIS(int t){
        
        if(t==0)
            return 0;
        else if(AIS_T==1 && t==1)
            return 1;
        
        double result=pow((t+0.0)/AIS_T, AIS_power);
        
        return result;
    }
    
    int MoveWithinAIS::moveTimes(int t, Param* param, RecTreeAux* rectreeAux, double* ll){
        
        RecTree* rectree=rectreeAux->rectree0;
        
        unsigned int start_0=rectreeAux->start, end_0=rectreeAux->end;
        unsigned int efrom_0=rectreeAux->efrom, eto_0=rectreeAux->eto;
        double tfromRel_0=rectreeAux->tfrom, ttoRel_0=rectreeAux->tto;
        double tfrom_0=tfromRel_0+rectree->getNode(efrom_0)->getAge(), tto_0=ttoRel_0+rectree->getNode(eto_0)->getAge();
        double l0=param->getLL();
        double TIMESTEPSIZE=0.1;
        double lpriorrat=0;// log of (prior ratio * transition ratio)
        
        vector<double> store_rectreeAux0=vector<double>(end_0-start_0);
        for (unsigned int i=start_0;i<end_0;i++)
            store_rectreeAux0[i-start_0]=rectreeAux->store_ll[i-start_0];
        
        RecEdge* edge_0 = new RecEdge(tfromRel_0,ttoRel_0,start_0,end_0,efrom_0,eto_0);
        lpriorrat-=rectree->priorEdge(edge_0,param);
        
        int movetto=gsl_rng_uniform_int(rng,2);//bernoulli: 1 moves tto, 0 moves tfrom
        
        int efrom=efrom_0;
        int eto=eto_0;
        double tfrom=tfrom_0;
        double tto=tto_0;
        double tmptime,tfromnew=tfrom,ttonew=tto;
        int start=start_0;
        int end=end_0;
        int tmpedge,etonew=eto,efromnew=efrom;
        double rootage=rectree->getNode(rectree->getN()*2-2)->getAge();
        
        dlog(1)<<"Proposing to change edge time "<<efrom<<":"<<tfrom<<"->"<<eto<<":"<<tto<<" to"<<flush;
        
        // change the edge
        if(movetto)
        {// moving the arrival time
            ttonew=tto+gsl_ran_gaussian(rng,TIMESTEPSIZE);
            while(ttonew<0 || ttonew>min(tfrom,rootage))
            {
                if(ttonew<0)
                    ttonew=-ttonew;
                if(ttonew>min(tfrom,rootage))
                    ttonew=2.0*min(tfrom,rootage)-ttonew;
            }// while loop as can bounce off reflecting boundaries several times
            tmptime=ttonew;
            tmpedge=etonew;
        }
        else
        {// moving the departure time
            tfromnew=tfrom+gsl_ran_gaussian(rng,TIMESTEPSIZE);
            if(tfromnew<tto)
                tfromnew=2.0*tto-tfromnew;
            tmptime=tfromnew;
            tmpedge=efromnew;
        }
        // update the edge index
        while(tmptime<rectree->getNode(tmpedge)->getAge())
        {
            if(gsl_rng_uniform(rng)<0.5)
                tmpedge=rectree->getNode(tmpedge)->getLeft()->getId();
            else
                tmpedge=rectree->getNode(tmpedge)->getRight()->getId();
            lpriorrat+=log(2.0);
        }
        if(tmpedge!=rectree->getN()*2-2)
        {
            while(tmptime>rectree->getNode(tmpedge)->getFather()->getAge())
            {
                tmpedge=rectree->getNode(tmpedge)->getFather()->getId();
                lpriorrat-=log(2.0);
                if(tmpedge==rectree->getN()*2-2)
                    break;
            }
        }
        if(movetto)
        {
            ttonew=tmptime;
            etonew=tmpedge;
        }
        else
        {
            tfromnew=tmptime;
            efromnew=tmpedge;
        }
        dlog(1)<<" "<<efromnew<<":"<<tfromnew<<"->"<<etonew<<":"<<ttonew<<"..."<<flush;
        
        double tfromnewRel=tfromnew-rectree->getNode(efromnew)->getAge(), ttonewRel=ttonew-rectree->getNode(etonew)->getAge();
        
        RecEdge* edge_new = new RecEdge(tfromnewRel,ttonewRel,start_0,end_0,efromnew,etonew);
        lpriorrat+=rectree->priorEdge(edge_new,param);
        
        rectreeAux->setAll(start,end,tfromnewRel,ttonewRel,efromnew,etonew);
        
        double ll_partial=rectreeAux->computePartialLL(param->getTheta());
        double ll0_partial=0;
        
        double l_star=l0;
        
        for (unsigned int i=start;i<end;i++)
            ll0_partial+=param->getLLsite(i);
        
        l_star-=ll0_partial;
        l_star+=ll_partial;
        
        delete edge_0;
        delete edge_new;
        // acceptance step
        
        if (log(gsl_rng_uniform(rng))>(l_star-(*ll))*gammaAIS(t)+lpriorrat)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            rectreeAux->setAll(start_0,end_0,tfromRel_0,ttoRel_0,efrom_0,eto_0);
            
            for (int i=start_0;i<end_0;i++)
                rectreeAux->store_ll[i-start_0]=store_rectreeAux0[i-start_0];
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        *ll=l_star;
        
        return(1);
    }
    
    double log_dbetaBin(int k, int n, double alpha, double beta){
        
        double result=gsl_sf_lngamma(n+1.0)-gsl_sf_lngamma(n-k+1.0) -gsl_sf_lngamma(k+1.0);
        
        result+=gsl_sf_lnbeta(k+alpha, n-k+beta);
        result-=gsl_sf_lnbeta(alpha, beta);
        
        return result;
    }
    
    int rbetaBin(int n, double alpha, double beta){
        
        double p=gsl_ran_beta(rng, alpha, beta);
        int result=gsl_ran_binomial(rng,p,n);
        
        return result;
    }
    
    int MoveWithinAIS::moveStartEnd(int t, Param* param, RecTreeAux* rectreeAux, double* ll){
        
        RecTree* rectree=rectreeAux->rectree0;
        
        unsigned int start_0=rectreeAux->start, end_0=rectreeAux->end;
        unsigned int efrom_0=rectreeAux->efrom, eto_0=rectreeAux->eto;
        double tfromRel_0=rectreeAux->tfrom, ttoRel_0=rectreeAux->tto;
        double l0=param->getLL();
        double lpriorrat=0;// log of (prior ratio * transition ratio)
        
        vector<double> store_rectreeAux0=vector<double>(end_0-start_0);
        for (unsigned int i=start_0;i<end_0;i++)
            store_rectreeAux0[i-start_0]=rectreeAux->store_ll[i-start_0];
        
        RecEdge* edge_0 = new RecEdge(tfromRel_0,ttoRel_0,start_0,end_0,efrom_0,eto_0);
        lpriorrat-=rectree->priorEdge(edge_0,param);
        
        int L=rectree->getL();
        int n=L-(end_0-start_0);
        double m=(start_0+1.0)/(L+1);
        double s=0.005*L;
        double alpha=s*m, beta=s*(1-m);
        
        int start_star=rbetaBin(n,alpha,beta);
        int end_star=end_0+start_star-start_0;
        
        lpriorrat-=log_dbetaBin(start_star, n, alpha, beta);
        
        m=(start_star+1.0)/(L+1);
        alpha=s*m;
        beta=s*(1-m);
        
        lpriorrat+=log_dbetaBin(start_0, n, alpha, beta);
        
        RecEdge* edge_star = new RecEdge(tfromRel_0,ttoRel_0,start_star,end_star,efrom_0,eto_0);
        lpriorrat+=rectree->priorEdge(edge_star,param);
        
        rectreeAux->setAll(start_star,end_star,tfromRel_0,ttoRel_0,efrom_0,eto_0);
        
        double ll_partial=rectreeAux->computePartialLL(param->getTheta());
        double ll0_partial=0;
        
        double l_star=l0;
        
        for (unsigned int i=start_star;i<end_star;i++)
            ll0_partial+=param->getLLsite(i);
        
        l_star-=ll0_partial;
        l_star+=ll_partial;
        
        delete edge_0;
        delete edge_star;
        // acceptance step
        
        if (log(gsl_rng_uniform(rng))>(l_star-(*ll))*gammaAIS(t)+lpriorrat)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            rectreeAux->setAll(start_0,end_0,tfromRel_0,ttoRel_0,efrom_0,eto_0);
            
            for (int i=start_0;i<end_0;i++)
                rectreeAux->store_ll[i-start_0]=store_rectreeAux0[i-start_0];
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        *ll=l_star;
        
        return(1);
    }
    
    int MoveWithinAIS::movePrior(int t, Param* param, RecTreeAux* rectreeAux, double* ll){
        
        RecTree* rectree=rectreeAux->rectree0;
        
        unsigned int start_0=rectreeAux->start, end_0=rectreeAux->end;
        unsigned int efrom_0=rectreeAux->efrom, eto_0=rectreeAux->eto;
        double tfromRel_0=rectreeAux->tfrom, ttoRel_0=rectreeAux->tto;
        double l0=param->getLL();
        
        vector<double> store_rectreeAux0=vector<double>(end_0-start_0);
        for (unsigned int i=start_0;i<end_0;i++)
            store_rectreeAux0[i-start_0]=rectreeAux->store_ll[i-start_0];
        
        
        double tfrom_star,tto_star;
        unsigned int start_star, end_star;
        unsigned int efrom_star, eto_star;
        double l_star;
        
        //Draw start and end
        rectree->setBlock(&start_star,&end_star,param->getDelta(),param->getData()->getBlocks());
        
        //Draw eto and tto
        eto_star=rectree->getPoint(&tto_star);
        
        //Draw efrom and tfrom
        tfrom_star=tto_star+rectree->getNode(eto_star)->getAge();
        efrom_star=rectree->getEdgeCoal(&tfrom_star);
        tfrom_star-=rectree->getNode(efrom_star)->getAge();
        
        rectreeAux->setAll(start_star,end_star,tfrom_star,tto_star,efrom_star,eto_star);
        
        double ll_partial=rectreeAux->computePartialLL(param->getTheta());
        double ll0_partial=0;
        for (unsigned int i=start_star;i<end_star;i++){
            ll0_partial+=param->getLLsite(i);
        }
        
        l_star=l0-ll0_partial+ll_partial;
        
        if (log(gsl_rng_uniform(rng))>(l_star-(*ll))*gammaAIS(t))
        {
            dlog(1)<<"Rejected!"<<endl;
            
            rectreeAux->setAll(start_0,end_0,tfromRel_0,ttoRel_0,efrom_0,eto_0);
            
            for (int i=start_0;i<end_0;i++)
                rectreeAux->store_ll[i-start_0]=store_rectreeAux0[i-start_0];
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        *ll=l_star;
        
        return(1);
    }
    
    int MoveWithinAIS::moveFull(int t, Param* param, RecTreeAux* rectreeAux, double* ll){
        
        RecTree* rectree=rectreeAux->rectree0;
        
        unsigned int start_0=rectreeAux->start, end_0=rectreeAux->end;
        unsigned int efrom_0=rectreeAux->efrom, eto_0=rectreeAux->eto;
        double tfromRel_0=rectreeAux->tfrom, ttoRel_0=rectreeAux->tto;
        double tfrom_0=tfromRel_0+rectree->getNode(efrom_0)->getAge(), tto_0=ttoRel_0+rectree->getNode(eto_0)->getAge();
        double l0=param->getLL();
        double TIMESTEPSIZE=0.1;
        double lpriorrat=0;// log of (prior ratio * transition ratio)
        
        vector<double> store_rectreeAux0=vector<double>(end_0-start_0);
        for (unsigned int i=start_0;i<end_0;i++)
            store_rectreeAux0[i-start_0]=rectreeAux->store_ll[i-start_0];
        
        RecEdge* edge_0 = new RecEdge(tfromRel_0,ttoRel_0,start_0,end_0,efrom_0,eto_0);
        lpriorrat-=rectree->priorEdge(edge_0,param);
        
        int movetto=gsl_rng_uniform_int(rng,2);//bernoulli: 1 moves tto, 0 moves tfrom
        
        int efrom=efrom_0;
        int eto=eto_0;
        double tfrom=tfrom_0;
        double tto=tto_0;
        double tmptime,tfromnew=tfrom,ttonew=tto;
        int tmpedge,etonew=eto,efromnew=efrom;
        double rootage=rectree->getNode(rectree->getN()*2-2)->getAge();
        
        dlog(1)<<"Proposing to change edge time "<<efrom<<":"<<tfrom<<"->"<<eto<<":"<<tto<<" to"<<flush;
        
        // change the edge
        if(movetto)
        {// moving the arrival time
            ttonew=tto+gsl_ran_gaussian(rng,TIMESTEPSIZE);
            while(ttonew<0 || ttonew>min(tfrom,rootage))
            {
                if(ttonew<0)
                    ttonew=-ttonew;
                if(ttonew>min(tfrom,rootage))
                    ttonew=2.0*min(tfrom,rootage)-ttonew;
            }// while loop as can bounce off reflecting boundaries several times
            tmptime=ttonew;
            tmpedge=etonew;
        }
        else
        {// moving the departure time
            tfromnew=tfrom+gsl_ran_gaussian(rng,TIMESTEPSIZE);
            if(tfromnew<tto)
                tfromnew=2.0*tto-tfromnew;
            tmptime=tfromnew;
            tmpedge=efromnew;
        }
        // update the edge index
        while(tmptime<rectree->getNode(tmpedge)->getAge())
        {
            if(gsl_rng_uniform(rng)<0.5)
                tmpedge=rectree->getNode(tmpedge)->getLeft()->getId();
            else
                tmpedge=rectree->getNode(tmpedge)->getRight()->getId();
            lpriorrat+=log(2.0);
        }
        if(tmpedge!=rectree->getN()*2-2)
        {
            while(tmptime>rectree->getNode(tmpedge)->getFather()->getAge())
            {
                tmpedge=rectree->getNode(tmpedge)->getFather()->getId();
                lpriorrat-=log(2.0);
                if(tmpedge==rectree->getN()*2-2)
                    break;
            }
        }
        if(movetto)
        {
            ttonew=tmptime;
            etonew=tmpedge;
        }
        else
        {
            tfromnew=tmptime;
            efromnew=tmpedge;
        }
        dlog(1)<<" "<<efromnew<<":"<<tfromnew<<"->"<<etonew<<":"<<ttonew<<"..."<<flush;
        
        double tfromnewRel=tfromnew-rectree->getNode(efromnew)->getAge(), ttonewRel=ttonew-rectree->getNode(etonew)->getAge();
        
        int L=rectree->getL();
        int n=L-(end_0-start_0);
        double m=(start_0+1.0)/(L+1);
        double s=0.005*L;
        double alpha=s*m, beta=s*(1-m);
        
        int start_star=rbetaBin(n,alpha,beta);
        int end_star=end_0+start_star-start_0;
        
        lpriorrat-=log_dbetaBin(start_star, n, alpha, beta);
        
        m=(start_star+1.0)/(L+1);
        alpha=s*m;
        beta=s*(1-m);
        
        lpriorrat+=log_dbetaBin(start_0, n, alpha, beta);
        
        RecEdge* edge_new = new RecEdge(tfromnewRel,ttonewRel,start_star,end_star,efromnew,etonew);
        lpriorrat+=rectree->priorEdge(edge_new,param);
        
        rectreeAux->setAll(start_star,end_star,tfromnewRel,ttonewRel,efromnew,etonew);
        
        double ll_partial=rectreeAux->computePartialLL(param->getTheta());
        double ll0_partial=0;
        
        double l_star=l0;
        
        for (unsigned int i=start_star;i<end_star;i++)
            ll0_partial+=param->getLLsite(i);
        
        l_star-=ll0_partial;
        l_star+=ll_partial;
        
        delete edge_0;
        delete edge_new;
        // acceptance step
        
        if (log(gsl_rng_uniform(rng))>(l_star-(*ll))*gammaAIS(t)+lpriorrat)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            rectreeAux->setAll(start_0,end_0,tfromRel_0,ttoRel_0,efrom_0,eto_0);
            
            for (int i=start_0;i<end_0;i++)
                rectreeAux->store_ll[i-start_0]=store_rectreeAux0[i-start_0];
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        *ll=l_star;
        
        return(1);
    }
    
    int MoveWithinAIS::moveFwd(int t, Param* param, RecTreeAux* rectreeAux, double* ll){
        
        int accept=0;
        
        //accept+=moveTimes(t,param,rectreeAux,ll);
        //accept+=moveStartEnd(t,param,rectreeAux,ll);
        //accept+=moveFull(t,param,rectreeAux,ll);
        accept+=movePrior(t,param,rectreeAux,ll);
        
        return accept;
    }
    
}
