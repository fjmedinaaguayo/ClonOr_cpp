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
        description= "Update removing a recombinant edge using AISRJ";
        name= "RemEdge";
        
    }
    
    double MoveRemEdgeAIS::gammaAIS(int t, int up){
        
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
    
    int MoveRemEdgeAIS::move(vector<int> * samplespace)
    {
        
        int T=param->getT_AIS();
        int N=param->getN_MAIS()+1;
        RecTree* rectree0=param->getRecTree();
        
        int which0=rectree0->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
        if (which0<0)
        {
            dlog(1)<<"No valid edges to change!"<<endl;
            return(-1);
        }
        
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
        
        double l0=param->getLL(), lu=log(gsl_rng_uniform(rng));
        
        for(int n=0; n<N; n++){
            
            if(n==0){
                
                for(int t=0; t<T; t++){
                    
                    if(t==0)
                        lratio[0]=l0-(gammaAIS(t+1,0)-gammaAIS(t,0))*ll[n]+log((1.0+rectree0->numRecEdge())*2.0/param->getRho()/rectree0->getTTotal());
                    else
                        lratio[n]+=-(gammaAIS(t+1,0)-gammaAIS(t,0))*ll[n];
                    
                    if(lu<=lratio[n]){
                        break;
                    }
                    
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
                        
                        if (log(gsl_rng_uniform(rng))<=(ll_Temp-ll[n])*gammaAIS(t+1,0)){
                            
                            dlog(1)<<"AISRJ t="<<t<<" out of "<<T<<"..."<<"MCMC in AIS Accepted!"<<endl;
                            ll[n]=ll_Temp;
                        }
                    }
                }
                
            }
        }
        
        numcalls++;
        dlog(1)<<"Proposing to remove edge via AISRJ "<<efrom[0]<<":"<<tfrom[0]<<"->"<<eto[0]<<":"<<tto[0]<<"...";
        
        if (lu>lratio[0])
        {
            dlog(1)<<"Rejected!"<<endl;
            
            which[0]=rectree0->addRecEdge_FMA(tfrom[0],tto[0],start[0],end[0],efrom[0],eto[0]);
            
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

    /*
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
                lratio+=l1-(gammaAIS(t+1)-gammaAIS(t))*l0+log((1.0+rectree->numRecEdge())*2.0/param->getRho()/rectree->getTTotal());
            else
                lratio+=-(gammaAIS(t+1)-gammaAIS(t))*l0;
            
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
                
                if (log(gsl_rng_uniform(rng))<=(l_star-l0)*(1-gammaAIS(t+1))){
                    
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
     */
    
    MoveRemEdgeAIS::~MoveRemEdgeAIS()
    {}
    //
    
} // end namespace weakarg
