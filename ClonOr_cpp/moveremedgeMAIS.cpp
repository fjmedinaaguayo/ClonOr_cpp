//
//  moveremedgeMAIS.cpp
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 30/04/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "moveremedgeMAIS.h"
#include <omp.h>
#include <algorithm>
//#include "rectreeAux.h"
//
using namespace std;
namespace weakarg
{
    
    MoveRemEdgeMAIS::MoveRemEdgeMAIS(Param *p,double a)
    : Move(p,a)
    {
        description= "Update removing a recombinant edge using MAISRJ";
        name= "RemEdge";
    }
    
    double MoveRemEdgeMAIS::gammaAIS(int t, int up){
        
        if(t==0)
            return 0;
        else if(param->getT_AIS()==1 && t==1)
            return 1;
        
        double result;
        if(up==1)
            result=pow((t+0.0)/param->getT_AIS() ,param->getgamma_AIS());
        else
            result=1-pow(1-(t+0.0)/param->getT_AIS() ,param->getgamma_AIS());
        
        return result;
    }
    
    double MoveRemEdgeMAIS::logSumExp(vector<double> x){
        
        double result=0;
        int n=x.size();
        
        double xmax=*max_element(x.begin(),x.end());
        
        for(int i=0; i<n; i++){
            x[i]-=xmax;
            result+=exp(x[i]);
        }
        result=log(result)+xmax;
        
        return(result);
    }
    
    vector<int> MoveRemEdgeMAIS::syst_resamp(vector<double> lw, int N){
        
        int sizeW=lw.size();
        vector<int> result(N);
        vector<double> lw_norm(sizeW);
        double lnormC=logSumExp(lw);
        
        for(int i=0; i<sizeW; i++){
            
            lw_norm[i]=lw[i]-lnormC;
        }
        
        double u=gsl_rng_uniform(rng);
        double C=0;
        int j=1;
        
        for (int i=0; i<N; i++){
            
            C+=exp(lw_norm[i]);
            while((u+j-1.0)/N<=C){
                
                result[j-1]=i;
                j+=1;
            }
        }
        
        return(result);
    }
    
    int MoveRemEdgeMAIS::mult_resamp(vector<double> lw){
        
        int sizeW=lw.size();
        int result;
        vector<double> lw_norm(sizeW);
        double lnormC=logSumExp(lw);
        
        for(int i=0; i<sizeW; i++){
            
            lw_norm[i]=lw[i]-lnormC;
        }
        
        double u=gsl_rng_uniform(rng);
        double C=0;
        int i=-1;
        
        while(u>C){
            
            i++;
            C+=exp(lw_norm[i]);
        }
        
        result=i;
        
        return(result);
    }
    
    int MoveRemEdgeMAIS::moveWithinAIS(int t, int up, RecTreeAux* rectreeAux, vector<double> store_ll, double* ll){
        
        RecTree* rectree=param->getRecTree();
        
        unsigned int start_0=rectreeAux->start, end_0=rectreeAux->end;
        unsigned int efrom_0=rectreeAux->efrom, eto_0=rectreeAux->eto;
        double tfromRel_0=rectreeAux->tfrom, ttoRel_0=rectreeAux->tto;
        double tfrom_0=tfromRel_0+rectree->getNode(efrom_0)->getAge(), tto_0=ttoRel_0+rectree->getNode(eto_0)->getAge();
        double l0=param->getLL();
        double TIMESTEPSIZE=0.1;
        double lpriorrat=0;// log of (prior ratio * transition ratio)
        
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
        vector<double> store(end-start);
        for (int i=start;i<end;i++)
            store[i-start]=param->getLLsite(i);
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
        
        ll_partial=rectreeAux->computePartialLL(param->getTheta());
        for (unsigned int i=start;i<end;i++)
            ll0_partial+=param->getLLsite(i);
        
        l_star-=ll0_partial;
        l_star+=ll_partial;
        
        delete edge_0;
        delete edge_new;
        // acceptance step
        if (log(gsl_rng_uniform(rng))>(l_star-l0)*(1-gammaAIS(t,up))+lpriorrat)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            rectreeAux->setAll(start_0,end_0,tfromRel_0,ttoRel_0,efrom_0,eto_0);
            
            if(up==1){
                for (int i=start_0;i<end_0;i++)
                    rectreeAux->store_ll[i-start_0]=store_ll[i-start_0];
            }
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        *ll=l_star;
        
        numaccept++;
        return(1);
    }
    
    int MoveRemEdgeMAIS::move(vector<int> * samplespace)
    {
        
        int T=param->getT_AIS();
        int N=param->getN_MAIS();
        RecTree* rectree0=param->getRecTree();
        
        int which0=rectree0->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
        if (which0<0)
        {
            dlog(1)<<"No valid edges to change!"<<endl;
            return(-1);
        }
        
        vector<double> lerror(N,0.0);
        vector<double> tfrom(N),tto(N);
        vector<unsigned int> start(N),end(N);
        vector<unsigned int> efrom(N),eto(N);
        vector<int> which(N);
        vector<vector<double> > store_ll0(N);
        vector<vector<double> > store_ll(N);
        double l, ll_ini=param->getLL();
        vector<double> ll(N), lratio(N,0.0);
        
        tfrom[0]=param->getRecTree()->getRecEdge(which0)->getTimeFrom();
        tto[0]=param->getRecTree()->getRecEdge(which0)->getTimeTo();
        start[0]=param->getRecTree()->getRecEdge(which0)->getStart();
        end[0]=param->getRecTree()->getRecEdge(which0)->getEnd();
        efrom[0]=param->getRecTree()->getRecEdge(which0)->getEdgeFrom();
        eto[0]=param->getRecTree()->getRecEdge(which0)->getEdgeTo();
        ll[0]=ll_ini;
        store_ll0[0]=vector<double>(end[0]-start[0]);
        for (unsigned int i=start[0];i<end[0];i++){
            store_ll0[0][i-start[0]]=param->getLLsite(i);
        }
        
        rectree0->remRecEdge(which0);
        param->computeLikelihood(start[0],end[0]);
        
        double l0=param->getLL(), lu=log(gsl_rng_uniform(rng)), lratio_avg=0;
        
        //#pragma omp parallel num_threads(n_threads)
        {
            
            //int id = omp_get_thread_num();
            //int total = omp_get_num_threads();
            
            //for(int i=0; i<N/n_threads; i++){
            //int n=i+id*N/n_threads;
            
        }
        
        #pragma omp parallel for
        for(int n=0; n<N; n++){
            
            int nThreads=omp_get_num_threads();
            if(n==0){
                
                param->getRectreeAux_vec()[n]->setAll(start[0],end[0],tfrom[0],tto[0],efrom[0],eto[0]);
                
                for(int t=0; t<T; t++){
                    
                    if(t==0)
                        lratio[0]=l0-(gammaAIS(t+1,0)-gammaAIS(t,0))*ll[n]+log((1.0+rectree0->numRecEdge())*2.0/param->getRho()/rectree0->getTTotal());
                    else
                        lratio[n]+=-(gammaAIS(t+1,0)-gammaAIS(t,0))*ll[n];
                    
                    if(t!=T-1){

///*
                        int accepted=moveWithinAIS(t+1, 0, param->getRectreeAux_vec()[n], store_ll[n], &ll[n]);
                        
                        if(accepted==1){
                            
                            //store_ll[n]=rectreeAux->store_ll;
                            //ll[n]=ll_Temp;
                        }
//*/
/*
                        double tfrom_Temp,tto_Temp;
                        unsigned int start_Temp, end_Temp;
                        unsigned int efrom_Temp, eto_Temp;
                        vector<double> store_ll0_Temp;
                        vector<double> store_ll_Temp;
                        double ll_Temp;
                        
                        //Draw start and end
                        param->getRectreeAux_vec()[n]->rectree0->setBlock(&start_Temp,&end_Temp,param->getDelta(),param->getData()->getBlocks());
                        
                        //Draw eto and tto
                        eto_Temp=param->getRectreeAux_vec()[n]->rectree0->getPoint(&tto_Temp);
                        
                        //Draw efrom and tfrom
                        tfrom_Temp=tto_Temp+param->getRectreeAux_vec()[n]->rectree0->getNode(eto_Temp)->getAge();
                        efrom_Temp=param->getRectreeAux_vec()[n]->rectree0->getEdgeCoal(&tfrom_Temp);
                        tfrom_Temp-=param->getRectreeAux_vec()[n]->rectree0->getNode(efrom_Temp)->getAge();
                        
                        param->getRectreeAux_vec()[n]->setAll(start_Temp,end_Temp,tfrom_Temp,tto_Temp,efrom_Temp,eto_Temp);
                        store_ll0_Temp=vector<double>(end_Temp-start_Temp);
                        
                        double ll_partial_Temp=param->getRectreeAux_vec()[n]->computePartialLL(param->getTheta());
                        double ll0_partial_Temp=0;
                        for (unsigned int i=start_Temp;i<end_Temp;i++){
                            store_ll0_Temp[i-start_Temp]=param->getLLsite(i);
                            ll0_partial_Temp+=store_ll0_Temp[i-start_Temp];
                        }
                        store_ll_Temp=param->getRectreeAux_vec()[n]->store_ll;
                        
                        ll_Temp=l0-ll0_partial_Temp+ll_partial_Temp;
                        
                        if (log(gsl_rng_uniform(rng))<=(ll_Temp-ll[n])*(1-gammaAIS(t+1,0))){
                            
                            dlog(1)<<"AISRJ t="<<t<<" out of "<<T<<"..."<<"MCMC in AIS Accepted!"<<endl;
                            ll[n]=ll_Temp;
                        }
 */
                    }
                }
                
            }
            else{
                
                int t=0;
                //Draw start and end
                param->getRectreeAux_vec()[n]->rectree0->setBlock(&start[n],&end[n],param->getDelta(),param->getData()->getBlocks());
                
                //Draw eto and tto
                eto[n]=param->getRectreeAux_vec()[n]->rectree0->getPoint(&tto[n]);
                
                //Draw efrom and tfrom
                tfrom[n]=tto[n]+param->getRectreeAux_vec()[n]->rectree0->getNode(eto[n])->getAge();
                efrom[n]=param->getRectreeAux_vec()[n]->rectree0->getEdgeCoal(&tfrom[n]);
                tfrom[n]-=param->getRectreeAux_vec()[n]->rectree0->getNode(efrom[n])->getAge();
                
                param->getRectreeAux_vec()[n]->setAll(start[n],end[n],tfrom[n],tto[n],efrom[n],eto[n]);
                store_ll0[n]=vector<double>(end[n]-start[n]);
                
                double ll_partial=param->getRectreeAux_vec()[n]->computePartialLL(param->getTheta());
                double ll0_partial=0;
                for (unsigned int i=start[n];i<end[n];i++){
                    store_ll0[n][i-start[n]]=param->getLLsite(i);
                    ll0_partial+=store_ll0[n][i-start[n]];
                }
                store_ll[n]=param->getRectreeAux_vec()[n]->store_ll;
                
                ll[n]=l0-ll0_partial+ll_partial;
                
                while(t<T){
                    
                    if(t==0) lratio[n]+=(gammaAIS(t+1,1)-gammaAIS(t,1))*ll[n]-l0+log(param->getRho()*rectree0->getTTotal()/2.0/(rectree0->numRecEdge()+1));
                    else
                        lratio[n]+=(gammaAIS(t+1,1)-gammaAIS(t,1))*ll[n];
                    
                    if(t!=T-1){

///*
                        int accepted=moveWithinAIS(t+1, 1, param->getRectreeAux_vec()[n], store_ll[n], &ll[n]);
                        
                        if(accepted==1){
                            
                            store_ll[n]=param->getRectreeAux_vec()[n]->store_ll;
                            //ll[n]=ll_Temp;
                            
                            tfrom[n]=param->getRectreeAux_vec()[n]->tfrom;
                            tto[n]=param->getRectreeAux_vec()[n]->tto;
                            start[n]=param->getRectreeAux_vec()[n]->start;
                            end[n]=param->getRectreeAux_vec()[n]->end;
                            efrom[n]=param->getRectreeAux_vec()[n]->efrom;
                            eto[n]=param->getRectreeAux_vec()[n]->eto;
                            
                        }
//*/
/*
                        double tfrom_Temp,tto_Temp;
                        unsigned int start_Temp, end_Temp;
                        unsigned int efrom_Temp, eto_Temp;
                        vector<double> store_ll0_Temp;
                        vector<double> store_ll_Temp;
                        double ll_Temp;
                        
                        //Draw start and end
                        param->getRectreeAux_vec()[n]->rectree0->setBlock(&start_Temp,&end_Temp,param->getDelta(),param->getData()->getBlocks());
                        
                        //Draw eto and tto
                        eto_Temp=param->getRectreeAux_vec()[n]->rectree0->getPoint(&tto_Temp);
                        
                        //Draw efrom and tfrom
                        tfrom_Temp=tto_Temp+param->getRectreeAux_vec()[n]->rectree0->getNode(eto_Temp)->getAge();
                        efrom_Temp=param->getRectreeAux_vec()[n]->rectree0->getEdgeCoal(&tfrom_Temp);
                        tfrom_Temp-=param->getRectreeAux_vec()[n]->rectree0->getNode(efrom_Temp)->getAge();
                        
                        param->getRectreeAux_vec()[n]->setAll(start_Temp,end_Temp,tfrom_Temp,tto_Temp,efrom_Temp,eto_Temp);
                        store_ll0_Temp=vector<double>(end_Temp-start_Temp);
                        
                        double ll_partial_Temp=param->getRectreeAux_vec()[n]->computePartialLL(param->getTheta());
                        double ll0_partial_Temp=0;
                        for (unsigned int i=start_Temp;i<end_Temp;i++){
                            store_ll0_Temp[i-start_Temp]=param->getLLsite(i);
                            ll0_partial_Temp+=store_ll0_Temp[i-start_Temp];
                        }
                        store_ll_Temp=param->getRectreeAux_vec()[n]->store_ll;
                        
                        ll_Temp=l0-ll0_partial_Temp+ll_partial_Temp;
                        
                        if (log(gsl_rng_uniform(rng))<=(ll_Temp-ll[n])*gammaAIS(t+1,1)){
                            
                            dlog(1)<<"AISRJ t="<<t<<" out of "<<T<<"..."<<"MCMC in AIS Accepted!"<<endl;
                            tfrom[n]=tfrom_Temp;
                            tto[n]=tto_Temp;
                            start[n]=start_Temp;
                            end[n]=end_Temp;
                            efrom[n]=efrom_Temp;
                            eto[n]=eto_Temp;
                            ll[n]=ll_Temp;
                            
                            store_ll[n]=store_ll_Temp;
                            store_ll0[n]=store_ll0_Temp;
                        }
*/
                    }
                    t++;
                }
            }
        }
        
        lratio[0]=-lratio[0];
        
        if(logSumExp(lerror)!=log(N))
            return -1;
        
        lratio_avg=logSumExp(lratio)-log(N);
        
        numcalls++;
        dlog(1)<<"Proposing to remove edge via MAISRJ "<<efrom[0]<<":"<<tfrom[0]<<"->"<<eto[0]<<":"<<tto[0]<<"...";
        
        if (lu>-lratio_avg)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            which[0]=rectree0->addRecEdge_FMA(tfrom[0],tto[0],start[0],end[0],efrom[0],eto[0]);
            //param->computeLikelihood(start[0],end[0]);
            //l=param->getLL();
            
            for (unsigned int i=start[0];i<end[0];i++)
                param->setlocLL(i,store_ll0[0][i-start[0]]);
            param->setLL(ll_ini);
            l=ll_ini;
            
            if(fabs(l-ll_ini)>1e-6){
                
                cout<<"Error removing edge, diff in ll: "<<(l-ll[0])<<endl;
                
                for(unsigned int site=start[0]; site<end[0]; site++){
                    if(fabs(store_ll[0][site-start[0]]-param->getLLsite(site))>1e-6){
                        cout<<"site "<<site<<", ll: "<<store_ll[0][site-start[0]]<<" vs "<<param->getLLsite(site)<<endl;
                    }
                }

                return -1;
            }
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        numaccept++;
        
        return(1);
    }
    
    MoveRemEdgeMAIS::~MoveRemEdgeMAIS()
    {}
    
} // end namespace weakarg
