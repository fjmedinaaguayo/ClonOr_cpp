//
//  moveaddedgeMAIS.cpp
//  ClonOr
//
//  Created by Felipe Medina Aguayo on 01/03/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "moveaddedgeMAIS.h"
//
using namespace std;
namespace weakarg
{
    
    MoveAddEdgeMAIS::MoveAddEdgeMAIS(Param *p,double a)
    : Move(p,a)
    {
        description= "Update adding a recombinant edge using AISRJ";
        name= "AddEdge";
    }
    
    double MoveAddEdgeMAIS::gammaAIS(int t, int T_AIS){
        
        if(t==0)
            return 0;
        
        double result=pow((t+0.0)/T_AIS ,5);
        return(result);
    }
    
    double MoveAddEdgeMAIS::logSumExp(vector<double> x){
        
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
        int T_AIS=param->getT_AIS();
        int N=100;
        RecTree* rectree0=param->getRecTree();
        
        vector<RecTree*> rectree(N);
        for(int i=0;i<N;i++){
            
            rectree[i]=new RecTree(*rectree0);
        }
        
        vector<double> tfrom(N),tto(N);
        vector<unsigned int> start(N),end(N);
        vector<unsigned int> efrom(N),eto(N);
        vector<int> which(N);
        vector<vector<double>> store_ll0(N);
        vector<vector<double>> store_ll(N);
        vector<double> l(N), lratio(N,0.0);
        
        double tfrom_star,tto_star;
        unsigned int start_star,end_star;
        unsigned int efrom_star,eto_star;
        int which_star;
        double l0=param->getLL(), l_star, lu=log(gsl_rng_uniform(rng)), lratio_avg=0;
        
        for(int n=0; n<N; n++){
            
            //Draw start and end
            rectree[n]->setBlock(&start[n],&end[n],param->getDelta(),param->getData()->getBlocks());
            
            //Draw eto and tto
            bool insamplespace=false;
            while(!insamplespace){
                eto[n]=rectree[n]->getPoint(&tto[n]);
                //Draw efrom and tfrom
                tfrom[n]=tto[n]+rectree[n]->getNode(eto[n])->getAge();
                efrom[n]=rectree[n]->getEdgeCoal(&tfrom[n]);
                tfrom[n]-=rectree[n]->getNode(efrom[n])->getAge();
                if(samplespace==NULL) insamplespace=true;
                else if (samplespace->size()==0) insamplespace=true;
                else
                {
                    for(unsigned int i=0;i<samplespace->size();i++) if((int)eto[n]==samplespace->at(i)||(int)efrom[n]==samplespace->at(i))
                    {
                        insamplespace=true;break;
                    }
                }
            }
            
            store_ll0[n].resize(end[n]-start[n]);
            for (unsigned int i=start[n];i<end[n];i++)
                store_ll0[n][i-start[n]]=param->getLLsite(i);
            param->setRecTree(rectree[n]);
            which[n]=rectree[n]->addRecEdge(tfrom[n],tto[n],start[n],end[n],efrom[n],eto[n]);
            if(which[n]<0) return(-1);
            param->computeLikelihood(start[n],end[n]);
            
            for(int t=0; t<T_AIS; t++){
                
                l[n]=param->getLL();
                
                if(t==0)
                    lratio[n]+=(gammaAIS(t+1,T_AIS)-gammaAIS(t,T_AIS))*l[n]-l0+log(param->getRho()*rectree[n]->getTTotal()/2.0/rectree[n]->numRecEdge());
                else
                    lratio[n]+=(gammaAIS(t+1,T_AIS)-gammaAIS(t,T_AIS))*l[n];
                
                if(t!=T_AIS-1){
                
                    RecTree* rectree_star=new RecTree(*rectree0);
                    
                    //Draw start and end
                    rectree_star->setBlock(&start_star,&end_star,param->getDelta(),param->getData()->getBlocks());
                    //Draw eto and tto
                    bool insamplespace=false;
                    while(!insamplespace){
                        eto_star=rectree_star->getPoint(&tto_star);
                        //Draw efrom and tfrom
                        tfrom_star=tto_star+rectree_star->getNode(eto_star)->getAge();
                        efrom_star=rectree_star->getEdgeCoal(&tfrom_star);
                        tfrom_star-=rectree_star->getNode(efrom_star)->getAge();
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
                    
                    store_ll[n].resize(end[n]-start[n]);
                    for (unsigned int i=start[n];i<end[n];i++)
                        store_ll[n][i-start[n]]=param->getLLsite(i);
                    
                    param->setRecTree(rectree0);
                    for (unsigned int i=start[n];i<end[n];i++)
                        param->setlocLL(i,store_ll0[n][i-start[n]]);
                    param->setLL(l0);
                    
                    vector<double> store_star(end_star-start_star);
                    for (unsigned int i=start_star;i<end_star;i++)
                        store_star[i-start_star]=param->getLLsite(i);
                    param->setRecTree(rectree_star);
                    which_star=rectree_star->addRecEdge(tfrom_star,tto_star,start_star,end_star,efrom_star,eto_star);
                    if(which_star<0) return(-1);
                    param->computeLikelihood(start_star,end_star);
                    l_star=param->getLL();
                    
                    //MCMC move
                    if (log(gsl_rng_uniform(rng))>(l_star-l[n])*(gammaAIS(t+1,T_AIS))){
                        
                        param->setRecTree(rectree[n]);
                        for (unsigned int i=start_star;i<end_star;i++)
                            param->setlocLL(i,store_star[i-start_star]);
                        for (unsigned int i=start[n];i<end[n];i++)
                            param->setlocLL(i,store_ll[n][i-start[n]]);
                        param->setLL(l[n]);
                    }
                    else{
                        
                        dlog(1)<<"AISRJ t="<<t<<" out of "<<T_AIS<<"..."<<"MCMC in AIS Accepted!"<<endl;
                        
                        tfrom[n]=tfrom_star;
                        tto[n]=tto_star;
                        start[n]=start_star;
                        end[n]=end_star;
                        efrom[n]=efrom_star;
                        eto[n]=eto_star;
                        l[n]=l_star;
                        which[n]=which_star;
                        
                        delete rectree[n];
                        rectree[n]=new RecTree (*rectree_star);
                        store_ll0[n]=store_star;
                        
                        store_ll[n].resize(end[n]-start[n]);
                        for (unsigned int i=start[n];i<end[n];i++){
                            store_ll[n][i-start[n]]=param->getLLsite(i);
                            param->setlocLL(i,store_ll0[n][i-start[n]]);
                        }
                        param->setLL(l0);
                    }
                    delete rectree_star;
                }
            }
        }
        
        lratio_avg=logSumExp(lratio)-log(N);

        int k=mult_resamp(lratio);
        numcalls++;
        dlog(1)<<"Proposing to add edge via AISRJ "<<efrom[k]<<":"<<tfrom[k]<<"->"<<eto[k]<<":"<<tto[k]<<"...";
        
        
        if (lu>lratio_avg)
        {
            dlog(1)<<"Rejected!"<<endl;
            param->setRecTree(rectree0);
            for (unsigned int i=start[k];i<end[k];i++)
                param->setlocLL(i,store_ll0[k][i-start[k]]);
            param->setLL(l0);
            
            for(int i=0; i<N; i++)
                delete rectree[i];

            return(0);
            //param->computeLikelihood(start,end);
        }
        else dlog(1)<<"Accepted!"<<endl;
        param->setRecTree(rectree[k]);
        for (unsigned int i=start[k];i<end[k];i++)
            param->setlocLL(i,store_ll[k][i-start[k]]);
        param->computeLikelihood(start[k],end[k]);
        param->setLL(l[k]);
        numaccept++;
        delete rectree0;
        for(int i=0; i<N; i++){
            if(i!=k)
                delete rectree[i];
        }
        
        return(1);
    }
    
    MoveAddEdgeMAIS::~MoveAddEdgeMAIS()
    {}
    //
    
} // end namespace weakarg
