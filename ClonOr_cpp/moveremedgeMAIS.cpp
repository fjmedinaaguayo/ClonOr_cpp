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
        
        MoveWithinAIS moveAIS(T, param->getgamma_AIS());
        
        #pragma omp parallel for
        for(int n=0; n<N; n++){
            
            int nThreads=omp_get_num_threads();
            if(n==0){
                
                param->getRectreeAux_vec()[n]->setAll(start[0],end[0],tfrom[0],tto[0],efrom[0],eto[0]);
                
                for(int t=0; t<T; t++){
                    
                    if(t==0)
                        lratio[0]=l0+log((1.0+rectree0->numRecEdge())*2.0/param->getRho()/rectree0->getTTotal());
                        
                    lratio[n]-=(1-moveAIS.gammaAIS(T-t-1)-(1-moveAIS.gammaAIS(T-t)))*ll[n];
                    
                    if(t!=T-1){

                        int accep=moveAIS.moveFwd(T-t-1, param, param->getRectreeAux_vec()[n], &ll[n]);
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
                    
                    if(t==0)
                        lratio[n]=-l0+log(param->getRho()*rectree0->getTTotal()/2.0/(rectree0->numRecEdge()+1));
                        
                        lratio[n]+=(moveAIS.gammaAIS(t+1)-moveAIS.gammaAIS(t))*ll[n];
                    
                    if(t!=T-1){

                        int accep=moveAIS.moveFwd(t+1, param, param->getRectreeAux_vec()[n], &ll[n]);
                        
                        if(accep!=0){
                            
                            store_ll[n]=param->getRectreeAux_vec()[n]->store_ll;
                            
                            tfrom[n]=param->getRectreeAux_vec()[n]->tfrom;
                            tto[n]=param->getRectreeAux_vec()[n]->tto;
                            start[n]=param->getRectreeAux_vec()[n]->start;
                            end[n]=param->getRectreeAux_vec()[n]->end;
                            efrom[n]=param->getRectreeAux_vec()[n]->efrom;
                            eto[n]=param->getRectreeAux_vec()[n]->eto;
                            
                        }
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
            //param->computeLikelihood();
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
