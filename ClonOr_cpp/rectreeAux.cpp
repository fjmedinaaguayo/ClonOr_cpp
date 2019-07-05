//
//  rectreeAux.cpp
//  ClonOr_cpp
//
//  Created by Felipe Medina Aguayo on 11/06/2019.
//  Copyright © 2019 Felipe Medina Aguayo. All rights reserved.
//

#include "rectreeAux.h"

namespace weakarg
{

    RecTreeAux::RecTreeAux(Data* data, RecTree* rectree0) : data(data), rectree0(rectree0){
        
        tabNode=std::vector<int>(2*data->getN()-1);
        tabRec=std::vector<int>(MAXNODES);
        tabFather=std::vector<int>(MAXNODES,-1);
        age=std::vector<double>(MAXNODES,0.0);
        tabEdge=std::vector<int>(MAXNODES);
        
        tabSons=std::vector<std::vector<int> >(MAXNODES,vector<int>(2,-1));
        tabSonsDist=std::vector<std::vector<double> >(MAXNODES,vector<double>(2,-1));
    }
    
    RecTreeAux::RecTreeAux(RecTreeAux* Copy) : data(Copy->data), rectree0(Copy->rectree0){
        
        this->tabNode=Copy->tabNode;
        this->tabRec=Copy->tabRec;
        this->tabFather=Copy->tabFather;
        this->age=Copy->age;
        this->tabEdge=Copy->tabEdge;
        
        this->tabSons=Copy->tabSons;
        this->tabSonsDist=Copy->tabSonsDist;
    }
    
    void RecTreeAux::setAll(unsigned int start, unsigned int end, double tfrom, double tto, unsigned int efrom, unsigned int eto){
     
        this->start=start;
        this->end=end;
        this->tfrom=tfrom;
        this->tto=tto;
        this->efrom=efrom;
        this->eto=eto;
        this->store_ll=vector<double>(end-start);
    }
    
    RecTreeAux::~RecTreeAux(){
     
        this->data=NULL;
        this->rectree0=NULL;
    }
    
    void RecTreeAux::computeLocalTree(double theta, unsigned int site, std::vector<RecEdge*> edge){
    
        int n=data->getN();
        double thetaPerSite=theta/data->getL();
        
        /*
        std::vector<int> tabNode(2*n-1);///<Used only in makeLocalTree (made global for speed)
        std::vector<int> tabRec(maxnodes);///<Used only in makeLocalTree (made global for speed)
        std::vector<int> tabFather(maxnodes,-1);///<Used only in makeLocalTree (made global for speed)
        std::vector<double> age(maxnodes,0.0);///<Used only in makeLocalTree (made global for speed)
        std::vector<int> tabEdge(maxnodes);///<Used only in makeLocalTree (made global for speed)
        */
         
        //Sort nodes and departure points in increasing order of age
        int indTabEdge=0;
        int aff=rectree0->affecting(site)+1;
        for (int i=0;i<n;i++)
            tabNode[i]=i;
        unsigned int indNode=n;
        unsigned int indEdge=0;
        bool isNode=true;
        
        /*
        tabSons=std::vector<std::vector<int> >(maxnodes,vector<int>(2,-1));
        tabSonsDist=std::vector<std::vector<double> >(maxnodes,vector<double>(2,-1));
        */
         
        for (int j=0;j<n-1+aff;j++)
        {
            while (indEdge<edge.size() && edge[indEdge]->affectsSite(site)==false)
                indEdge++;
            if (indNode>=rectree0->nodes.size())
            {
                isNode=false;
                indEdge++;
                goto next;
            };
            if (indEdge>=edge.size())
            {
                isNode=true;
                indNode++;
                goto next;
            }
            if (rectree0->nodes[indNode]->getAge()<edge[indEdge]->getTimeFrom()+rectree0->nodes[edge[indEdge]->getEdgeFrom()]->getAge())
            {
                isNode=true;
                indNode++;
            }
            else
            {
                isNode=false;
                indEdge++;
            }
        next:
            if (isNode)
            {
                tabNode[indNode-1]=j+n;
                age[j+n]=rectree0->nodes[indNode-1]->getAge();
            }
            else
            {
                tabRec[indEdge-1]=j+n;
                age[j+n]=edge[indEdge-1]->getTimeFrom()+rectree0->nodes[edge[indEdge-1]->getEdgeFrom()]->getAge();
                //age[j+n]=rectree0->getEdgeTimeAbsFrom(indEdge-1);
                tabEdge[indTabEdge++]=indEdge-1;
            }
        }
        
        //Find children of each node
        for (int i=n;i<2*n-1+aff;i++)
            for (int isLeft=0;isLeft<=1;isLeft++)
            {
                //Find out which edge "ed" to go down from and starting from age "ageCurrent"
                int ind,ii;
                double ageCurrent;
                int ed;
                if (i<2*n-1)
                {
                    isNode=true;
                    ind=i;
                    ii=tabNode[ind];
                    ageCurrent=rectree0->nodes[ind]->getAge();
                    if (isLeft)
                    {
                        ed=rectree0->nodes[ind]->getLeft ()->getId();
                    }
                    else
                    {
                        ed=rectree0->nodes[ind]->getRight()->getId();
                    }
                }
                else
                {
                    isNode=false;
                    ind=tabEdge[i-2*n+1];
                    ii=tabRec[ind];
                    if (isLeft)
                    {
                        ed=edge[ind]->getEdgeFrom();
                        //ageCurrent=rectree0->getEdgeTimeAbsFrom(ind);
                        ageCurrent=edge[ind]->getTimeFrom()+rectree0->nodes[edge[ind]->getEdgeFrom()]->getAge();
                    }
                    else
                    {
                        ed=edge[ind]->getEdgeTo();
                        //ageCurrent=rectree0->getEdgeTimeAbsTo(ind);
                        ageCurrent=edge[ind]->getTimeTo()+rectree0->nodes[edge[ind]->getEdgeTo()]->getAge();
                    }
                }
                
                //Go down the edge to find the son (if any)
                int cur=-1;
                double curAge=0;
                bool from=false;
                for (int k=0;k<aff;k++)
                {
                    if (edge[tabEdge[k]]->getEdgeFrom()==ed && isless(edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge(),ageCurrent) && (cur==-1||isless(curAge,edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge())))
                    {
                        cur=tabEdge[k];
                        //curAge=rectree0->getEdgeTimeAbsFrom(tabEdge[k]);
                        curAge=edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge();
                        from=true;
                    };
                    if (edge[tabEdge[k]]->getEdgeTo()==ed && isless(edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge(),ageCurrent) && (cur==-1||isless(curAge,edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge())))
                    {
                        cur=tabEdge[k];
                        //curAge=rectree0->getEdgeTimeAbsTo(tabEdge[k]);
                        curAge=edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge();
                        from=false;
                    };
                }
                if (cur==-1)
                {
                    isNode=true;
                    ind=ed;
                }
                else if (from)
                {
                    isNode=false;
                    ind=cur;
                }
                else
                {
                    tabSons[ii][isLeft]=-1;
                    continue;
                }
                if (isNode)
                    tabSons[ii][isLeft]=tabNode[ind];
                else
                    tabSons[ii][isLeft]=tabRec[ind];
                tabFather[tabSons[ii][isLeft]]=ii;
            }
        
        //Remove unnecessary nodes
        int cur=n;
        for (int i=n;i<2*n-1+aff;i++)
        {
            int f=tabFather[i];
            if (tabSons[i][0]==-1 && tabSons[i][1]==-1)
            {
                if (tabSons[f][0]==i)
                    tabSons[f][0]=-1;
                else
                    tabSons[f][1]=-1;
                continue;
            }
            if (tabSons[i][0]==-1 && f>=0)
            {
                if (tabSons[f][0]==i)
                    tabSons[f][0]=tabSons[i][1];
                else
                    tabSons[f][1]=tabSons[i][1];
                tabSons[i][1]=-1;
                continue;
            };
            if (tabSons[i][1]==-1 && f>=0)
            {
                if (tabSons[f][0]==i)
                    tabSons[f][0]=tabSons[i][0];
                else
                    tabSons[f][1]=tabSons[i][0];
                tabSons[i][0]=-1;
                continue;
            };
            if (tabSons[i][0]>=0 && tabSons[i][1]>=0)
            {
                tabSons[cur][0]=tabSons[i][0];
                tabSons[cur][1]=tabSons[i][1];
                age[cur]=age[i];
                if (f>=0)
                {
                    if (tabSons[f][0]==i)
                        tabSons[f][0]=cur;
                    else
                        tabSons[f][1]=cur;
                }
                
                tabSonsDist[cur][0]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][0]])*thetaPerSite);
                tabSonsDist[cur][1]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][1]])*thetaPerSite);
                if(tabSonsDist[cur][0]>1||tabSonsDist[cur][1]>1)
                {
                    cerr<<std::setprecision(64)<<"makeLocalTree:negative sonsdist! "<<age[cur]<<"-"<<age[tabSons[cur][0]]<<" and "<<age[cur]<<"-"<<age[tabSons[cur][1]]<<endl;
                    throw("Negative sons dist");
                }
                cur++;
            }
        }
        
    }
    
    double RecTreeAux::computePartialLL(double theta){
    
        int n=data->getN();
        double ll_Partial=0;
        
        vector<vector<double> > f=vector<vector<double> >(n*2-1,vector<double>(4,0.0)); //vector for likelihood
        
        unsigned int i=0;
        while (i<rectree0->edge.size() && rectree0->edge[i]->getTimeFrom()+rectree0->nodes[rectree0->edge[i]->getEdgeFrom()]->getAge()<tfrom+rectree0->nodes[efrom]->getAge())
            i++;
        
        RecEdge* re = new RecEdge(tfrom,tto,start,end,efrom,eto);
        
        std::vector<RecEdge*> edge(rectree0->edge.size());
        for(int i=0; i<edge.size(); i++){
            
            edge[i]=new RecEdge(rectree0->edge[i]->timeFrom,rectree0->edge[i]->timeTo,rectree0->edge[i]->gstart,rectree0->edge[i]->gend,rectree0->edge[i]->edgeFrom,rectree0->edge[i]->edgeTo);
        }
        edge.insert(edge.begin()+i,re);
         
        for (unsigned int site=start;site<end;site++){
            
            //make local tree at start and when other recombination event affects
            if(site==start || !rectree0->sameLocalTreeAsPrev(site)){
            //if(site>=start){
                
                computeLocalTree(theta,site,edge);
                
            /*
                //Sort nodes and departure points in increasing order of age
                int indTabEdge=0;
                int aff=rectree0->affecting(site)+1;
                for (int i=0;i<n;i++)
                    tabNode[i]=i;
                unsigned int indNode=n;
                unsigned int indEdge=0;
                bool isNode=true;
                
                for (int j=0;j<n-1+aff;j++)
                {
                    while (indEdge<edge.size() && edge[indEdge]->affectsSite(site)==false)
                        indEdge++;
                    if (indNode>=rectree0->nodes.size())
                    {
                        isNode=false;
                        indEdge++;
                        goto next;
                    };
                    if (indEdge>=edge.size())
                    {
                        isNode=true;
                        indNode++;
                        goto next;
                    }
                    if (rectree0->nodes[indNode]->getAge()<edge[indEdge]->getTimeFrom()+rectree0->nodes[edge[indEdge]->getEdgeFrom()]->getAge())
                    {
                        isNode=true;
                        indNode++;
                    }
                    else
                    {
                        isNode=false;
                        indEdge++;
                    }
                next:
                    if (isNode)
                    {
                        tabNode[indNode-1]=j+n;
                        age[j+n]=rectree0->nodes[indNode-1]->getAge();
                    }
                    else
                    {
                        tabRec[indEdge-1]=j+n;
                        age[j+n]=edge[indEdge-1]->getTimeFrom()+rectree0->nodes[edge[indEdge-1]->getEdgeFrom()]->getAge();
                        //age[j+n]=rectree0->getEdgeTimeAbsFrom(indEdge-1);
                        tabEdge[indTabEdge++]=indEdge-1;
                    }
                }
                
                //Find children of each node
                for (int i=n;i<2*n-1+aff;i++)
                    for (int isLeft=0;isLeft<=1;isLeft++)
                    {
                        //Find out which edge "ed" to go down from and starting from age "ageCurrent"
                        int ind,ii;
                        double ageCurrent;
                        int ed;
                        if (i<2*n-1)
                        {
                            isNode=true;
                            ind=i;
                            ii=tabNode[ind];
                            ageCurrent=rectree0->nodes[ind]->getAge();
                            if (isLeft)
                            {
                                ed=rectree0->nodes[ind]->getLeft ()->getId();
                            }
                            else
                            {
                                ed=rectree0->nodes[ind]->getRight()->getId();
                            }
                        }
                        else
                        {
                            isNode=false;
                            ind=tabEdge[i-2*n+1];
                            ii=tabRec[ind];
                            if (isLeft)
                            {
                                ed=edge[ind]->getEdgeFrom();
                                //ageCurrent=rectree0->getEdgeTimeAbsFrom(ind);
                                ageCurrent=edge[ind]->getTimeFrom()+rectree0->nodes[edge[ind]->getEdgeFrom()]->getAge();
                            }
                            else
                            {
                                ed=edge[ind]->getEdgeTo();
                                //ageCurrent=rectree0->getEdgeTimeAbsTo(ind);
                                ageCurrent=edge[ind]->getTimeTo()+rectree0->nodes[edge[ind]->getEdgeTo()]->getAge();
                            }
                        }
                        
                        //Go down the edge to find the son (if any)
                        int cur=-1;
                        double curAge=0;
                        bool from=false;
                        for (int k=0;k<aff;k++)
                        {
                            if (edge[tabEdge[k]]->getEdgeFrom()==ed && isless(edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge(),ageCurrent) && (cur==-1||isless(curAge,edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge())))
                            {
                                cur=tabEdge[k];
                                //curAge=rectree0->getEdgeTimeAbsFrom(tabEdge[k]);
                                curAge=edge[tabEdge[k]]->getTimeFrom()+rectree0->nodes[edge[tabEdge[k]]->getEdgeFrom()]->getAge();
                                from=true;
                            };
                            if (edge[tabEdge[k]]->getEdgeTo()==ed && isless(edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge(),ageCurrent) && (cur==-1||isless(curAge,edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge())))
                            {
                                cur=tabEdge[k];
                                //curAge=rectree0->getEdgeTimeAbsTo(tabEdge[k]);
                                curAge=edge[tabEdge[k]]->getTimeTo()+rectree0->nodes[edge[tabEdge[k]]->getEdgeTo()]->getAge();
                                from=false;
                            };
                        }
                        if (cur==-1)
                        {
                            isNode=true;
                            ind=ed;
                        }
                        else if (from)
                        {
                            isNode=false;
                            ind=cur;
                        }
                        else
                        {
                            tabSons[ii][isLeft]=-1;
                            continue;
                        }
                        if (isNode)
                            tabSons[ii][isLeft]=tabNode[ind];
                        else
                            tabSons[ii][isLeft]=tabRec[ind];
                        tabFather[tabSons[ii][isLeft]]=ii;
                    }
                
                //Remove unnecessary nodes
                int cur=n;
                for (int i=n;i<2*n-1+aff;i++)
                {
                    int f=tabFather[i];
                    if (tabSons[i][0]==-1 && tabSons[i][1]==-1)
                    {
                        if (tabSons[f][0]==i)
                            tabSons[f][0]=-1;
                        else
                            tabSons[f][1]=-1;
                        continue;
                    }
                    if (tabSons[i][0]==-1 && f>=0)
                    {
                        if (tabSons[f][0]==i)
                            tabSons[f][0]=tabSons[i][1];
                        else
                            tabSons[f][1]=tabSons[i][1];
                        tabSons[i][1]=-1;
                        continue;
                    };
                    if (tabSons[i][1]==-1 && f>=0)
                    {
                        if (tabSons[f][0]==i)
                            tabSons[f][0]=tabSons[i][0];
                        else
                            tabSons[f][1]=tabSons[i][0];
                        tabSons[i][0]=-1;
                        continue;
                    };
                    if (tabSons[i][0]>=0 && tabSons[i][1]>=0)
                    {
                        tabSons[cur][0]=tabSons[i][0];
                        tabSons[cur][1]=tabSons[i][1];
                        age[cur]=age[i];
                        if (f>=0)
                        {
                            if (tabSons[f][0]==i)
                                tabSons[f][0]=cur;
                            else
                                tabSons[f][1]=cur;
                        }
                        
                        tabSonsDist[cur][0]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][0]])*thetaPerSite);
                        tabSonsDist[cur][1]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][1]])*thetaPerSite);
                        if(tabSonsDist[cur][0]>1||tabSonsDist[cur][1]>1)
                        {
                            cerr<<std::setprecision(64)<<"makeLocalTree:negative sonsdist! "<<age[cur]<<"-"<<age[tabSons[cur][0]]<<" and "<<age[cur]<<"-"<<age[tabSons[cur][1]]<<endl;
                            throw("Negative sons dist");
                        }
                        cur++;
                    }
                }
             */
            }
            
            //computeSiteLL
            
            for (int i=n-1;i>=0;i--)
            {
                if (data->get(i,site)>3)
                    for (int a=0;a<4;a++)
                        f[i][a]=1.0;
                else
                {
                    f[i][0]=0.0;
                    f[i][1]=0.0;
                    f[i][2]=0.0;
                    f[i][3]=0.0;
                    f[i][data->get(i,site)]=1.0;
                }
            }
            int y,z; // children of this site
            double ey,ez,eyeq,ezeq,eyneq,ezneq;
            for (int i=n;i<2*n-1;i++)
            {
                y=tabSons[i][true];
                ey=tabSonsDist[i][true];
                z=tabSons[i][false];
                ez=tabSonsDist[i][false];
                
                if(y<0||z<0||y>=(int)f.size()||z>=(int)f.size()){
                    // periodically warn the user about these problems, not every time it happens.
                    static long badTreeWarning=0;
                    if(badTreeWarning <= 100){
                        cerr<<"Error in computeSiteLL: Invalid Local Tree, happened " << badTreeWarning << " times" <<endl;
                        if(badTreeWarning==100)  cerr << "Not reporting any more errors, fix this program!!\n";
                    }
                    badTreeWarning++;
                    store_ll[site-start]=0;                 /// WARNING!  This does NOT throw an error but ignores the site!
                    
                }
                eyeq =0.25*(1.0+3.0*ey);
                ezeq =0.25*(1.0+3.0*ez);
                eyneq=0.25*(1.0-ey);
                ezneq=0.25*(1.0-ez);
                for (unsigned int j=0;j<4;j++)
                {
                    int j1=(j+1)%4,j2=(j+2)%4,j3=(j+3)%4;
                    double YnotJ=f[y][j1]+f[y][j2]+f[y][j3];
                    double ZnotJ=f[z][j1]+f[z][j2]+f[z][j3];
                    f[i][j]=(f[z][j]*ezeq+ZnotJ*ezneq)*(f[y][j]*eyeq+YnotJ*eyneq);
                }
            }
            store_ll[site-start]=LOG025+log(f[2*n-2][0]+f[2*n-2][1]+f[2*n-2][2]+f[2*n-2][3]);
            ll_Partial+=store_ll[site-start];
            
            if(isnan(store_ll[site-start]))
            {
                cerr<<"ERROR: computeSiteLL:Site "<<site<<"has NaN lik!"<<endl;
                for(int i=n;i<2*n-1;i++)
                {
                    cout<<"i="<<i<<endl;
                    for (unsigned int j=0;j<4;j++)
                        cout<<"j="<<j<<" f[i][j]="<<f[i][j]<<endl;
                }
                cout<<"locll["<<site<<"]=log(0.25)+log("<<f[2*n-2][0]+f[2*n-2][1]+f[2*n-2][2]+f[2*n-2][3]<<")"<<endl;
                throw("NaN lik");
            }
        }

        for(int i=0; i<edge.size(); i++){
            
            delete edge[i];
        }
        
        return ll_Partial;
    }
}
/*d
using namespace std;
namespace weakarg
{
    
    inline double fround(double n, double d)
    {
        return floor(n * pow(10., d) + .5) / pow(10., d);
    }
}

// setBlock assigns start stie and end site on sequence
// getPoint assigns time_to and edge_to using uniform distn on tree
// getNode helps to obtain node from edge
// getAge obtains age of node
// getEdgeCoal for obtaning edge from
// getNode used again for modifying tfrom



rectree_vec[n]->setBlock(&start[n],&end[n],param->getDelta(),param->getData()->getBlocks());
eto[n]=rectree_vec[n]->getPoint(&tto[n]);
tfrom[n]=tto[n]+rectree_vec[n]->getNode(eto[n])->getAge();
efrom[n]=rectree_vec[n]->getEdgeCoal(&tfrom[n]);
tfrom[n]-=rectree_vec[n]->getNode(efrom[n])->getAge();
which[n]=rectree_vec[n]->addRecEdge_FMA(tfrom[n],tto[n],start[n],end[n],efrom[n],eto[n]);
param_vec[n]->computeLikelihood(start[n],end[n]);
l[n]=param_vec[n]->getLL();

vector<double> temp(end[n]-start[n]);
for (unsigned int i=start[n];i<end[n];i++)
temp[i-start[n]]=param->getLLsite(i);

store_ll0[n]=temp;


lratio[n]+=(gammaAIS(1,T_AIS)-gammaAIS(0,T_AIS))*l[n]-l0+log(param_vec[n]->getRho()*rectree_vec[n]->getTTotal()/2.0/rectree_vec[n]->numRecEdge());
}
}

*/
