//
//  moveaddedgeMAIS.cpp
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 01/03/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "moveaddedgeMAIS.h"
#include <omp.h>
#include <algorithm>
//#include "rectreeAux.h"
//
using namespace std;
namespace weakarg
{
    
    MoveAddEdgeMAIS::MoveAddEdgeMAIS(Param *p,double a)
    : Move(p,a)
    {
        description= "Update adding a recombinant edge using MAISRJ";
        name= "AddEdge";
    }
    
    double MoveAddEdgeMAIS::gammaAIS(int t){
        
        if(t==0)
            return 0;
        
        double result=pow((t+0.0)/param->getT_AIS() , param->getgamma_AIS());
        return(result);
    }
    
    double MoveAddEdgeMAIS::logSumExp(vector<double> x){
        
        double result=0;
        int n=x.size();
        
        double xmax=*(std::max_element(x.begin(),x.end()));
        
        for(int i=0; i<n; i++){
            x[i]-=xmax;
            result+=exp(x[i]);
        }
        result=log(result)+xmax;
        
        return(result);
    }
    
    vector<int> MoveAddEdgeMAIS::syst_resamp(vector<double> lw, int N){
        
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
    
    int MoveAddEdgeMAIS::mult_resamp(vector<double> lw){
    
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
    
    int MoveAddEdgeMAIS::move(vector<int> * samplespace)
    {
        
        int T=param->getT_AIS();
        int N=param->getN_MAIS();
        RecTree* rectree0=param->getRecTree();
        
        vector<double> lerror(N,0.0);
        vector<double> tfrom(N),tto(N);
        vector<unsigned int> start(N),end(N);
        vector<unsigned int> efrom(N),eto(N);
        int which;
        vector<vector<double> > store_ll0(N);
        vector<vector<double> > store_ll(N);
        double l;
        vector<double> ll(N), lratio(N,0.0);
        
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
            int t=0;
            
            //Draw start and end
            param->getRectreeAux_vec()[n]->rectree0->setBlock(&start[n],&end[n],param->getDelta(),param->getData()->getBlocks());
            //rectree_vec[n]->setBlock(&start[n],&end[n],param->getDelta(),param->getData()->getBlocks());
            
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
                
                if(t==0) lratio[n]+=(gammaAIS(t+1)-gammaAIS(t))*ll[n]-l0+log(param->getRho()*rectree0->getTTotal()/2.0/(rectree0->numRecEdge()+1));
                else
                    lratio[n]+=(gammaAIS(t+1)-gammaAIS(t))*ll[n];
                
                if(t!=T-1){
                    
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
                    
                    if (log(gsl_rng_uniform(rng))<=(ll_Temp-ll[n])*(gammaAIS(t+1))){
                        
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
                    
                }
                t++;
                
            }
        }
        
        if(logSumExp(lerror)!=log(N))
            return -1;
        
        lratio_avg=logSumExp(lratio)-log(N);

        int k=mult_resamp(lratio);
        numcalls++;
        dlog(1)<<"Proposing to add edge via MAISRJ "<<efrom[k]<<":"<<tfrom[k]<<"->"<<eto[k]<<":"<<tto[k]<<"...";
        
        
        if (lu>lratio_avg)
        {
            dlog(1)<<"Rejected!"<<endl;
            
            return(0);
        }
        else dlog(1)<<"Accepted!"<<endl;
        
        which=rectree0->addRecEdge_FMA(tfrom[k],tto[k],start[k],end[k],efrom[k],eto[k]);
        //param->computeLikelihood(start[k],end[k]);
        //l=param->getLL();
        
        for (unsigned int i=start[k];i<end[k];i++)
            param->setlocLL(i,store_ll[k][i-start[k]]);
        
        param->setLL(ll[k]);
        l=ll[k];
        
        if(fabs(l-ll[k])>1e-6){
            
            cout<<"Error adding edge, diff in ll: "<<(l-ll[k])<<endl;
            
            for(unsigned int site=start[k]; site<end[k]; site++){
                if(fabs(store_ll[k][site-start[k]]-param->getLLsite(site))>1e-6){
                    cout<<"site "<<site<<", ll: "<<store_ll[k][site-start[k]]<<" vs "<<param->getLLsite(site)<<endl;
                }
            }
            
            return -1;
        }
        numaccept++;
        
        return(1);
    }
    
    MoveAddEdgeMAIS::~MoveAddEdgeMAIS()
    {}
    //
    
} // end namespace weakarg
