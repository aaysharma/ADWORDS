#include<iostream>
#include<vector>
#include<iterator>
#include<random>
#include<cmath>
#include<algorithm>
#include<climits>
using namespace std;


class OnlineBipartite{
public:
	int m, n;						// two parititions adverisers A and Queries Q, |A|=m, |Q|=n
	vector<vector<double>> bids;		// adjacency list of bids b/w A-Q
	vector<double> B;					// budgets of A Bj
	vector<double> w;                  // uniformly randomly assigned weights to A wj
	vector<double> p;                  // prices of A pj
	vector<double> r;                  // revenue of A rj
	vector<double> L;                  // limit of A Lj

	OnlineBipartite(int a, int q) {
		m=a; n=q;
        bids.resize(m);				// resize the vector to hold `m` elements of type `vector<int>`
        for(int j=0; j<m; j++)
        	bids[j].resize(n);
        B.resize(m);
        w.resize(m);
        p.resize(m);
        r.resize(m);
        L.resize(m);
    }

    void addEdges(vector<double> edges, int i){// add edges for each new query q
        for (int j=0; j<m; j++)
            bids[j][i] = edges[j];
    }

    void addEdge(int j, double bid, int i){// add edge (a,q) with weightbid(a,q)
        bids[j][i] = bid;
    }

    void printGraph(){
        cout<<"Bids:"<<endl;
    	for(int j=0; j<m; j++){
    	    cout<<"bid("<<j<<", *): ";
    		for(int i=0; i<n; i++){
    		    if(i==n-1)
    		    	cout<<bids[j][i];
    		    else
    		    	cout<<bids[j][i]<<", ";
    		}
    		cout<<endl;
    	}
    	cout<<"Budgets[j]::"<<endl;
    	for(int j=0; j<m; j++)
    	    cout<<" "<<B[j];
    	cout<<endl;
    }
    void printAdj(){
        cout<<"Adjacency Matrix:"<<endl;
        
        for(int i=0; i<n+m; i++){
            for(int j=0; j<n+m; j++){
                if(i<n){
                    if(j<n){            // both in Q
                        cout<<0<<", ";
                        continue;
                    }
                    else{               // i in Q, j in A
                        cout<<bids[j%m][i]<<", ";
                    }
                }
                else{
                    if(j<n){        // i in A, j in Q
                        cout<<bids[i%m][j]<<", ";
                        continue;
                    }
                    else{           // i in A, j in A
                        cout<<0<<", ";
                    }
                }
            }
            cout<<endl;
        }

    	cout<<"Budget[j+m]::"<<endl;
    	for(int j=0; j<m; j++)
    	    cout<<B[j]<<" ";
    	cout<<endl;
    }
};

bool valid_matching(OnlineBipartite &G, vector<int> &M){
    int n=G.n, m=G.m;
    vector<double> cap(m,0);
    for(int i=0; i<n; i++){
        if(M[i]!=-1){
            cap[M[i]]+=G.bids[M[i]][i];
            if(cap[M[i]] > G.B[M[i]]) return 0;
        }
    }
    return 1;
}

double opt(OnlineBipartite &G, vector<int> &M, int curr){
    if(curr==G.n){    // after last index
        if(!valid_matching(G, M))
            return INT_MIN;
        return 0;
    }
    double mx = opt(G, M, curr+1);
    for(int j=0; j<G.m; j++){
        M[curr] = j;
        int x = G.bids[j][curr] + opt(G, M, curr+1);
        mx = (mx>x)?mx:x;
        // cout<<"curr= "<<curr<<" mx= "<<mx<<" j= "<<j<<endl;
    }
    
    return mx;
}


void test(){
    int U=20, L=2, U_bids=40;
    random_device rd;
	mt19937 gen0(rd());        // random uniform generator
	uniform_int_distribution<int> distr0(L, U); // sample m,n
	int m = distr0(gen0), n =distr0(gen0);
	cout<<"m="<<m<<" n="<<n<<endl;
	uniform_int_distribution<int> distr1(0, U_bids); // sample bids
	uniform_int_distribution<int> distr2(0, n*U_bids); // sample budgets


	OnlineBipartite G(m, n);

	for(int j=0; j<m; j++){					// read the input graph
		for(int i=0; i<n; i++){
			double bid = distr1(gen0);
			G.addEdge(j, bid, i);
		}
		double bj = distr2(gen0);
		G.B[j] = bj;
		
	}
	
	// main algorithm
	
	vector<int> M(n,-1);         // output matching
	double algo=0;
	double W=0, Wf=0;
	vector<int> u(n,0);
	default_random_engine generator;        // random uniform generator for Bj
	uniform_real_distribution<double> distr(0,1.0); // -do-
	
	
	for(int j=0; j<m; j++){					// update budgets of A and assign random prices
		G.w[j] = distr(generator);          // sample wj randomly
		G.p[j] = exp(G.w[j]-1);
		G.r[j] = 0;
		G.L[j] = G.B[j];
	}
	
	G.printGraph();
// 	G.printAdj();
	
	for(int i=0; i<n; i++){                 // queries in online input
	    vector<double> ebid(m);             // vector of effective bids
	    int argmax = -1;
	    for(int j=0; j<m; j++){
	       // if(G.L[j]>=G.bids[j][i]){
	        if(G.L[j]>0){
	            ebid[j] = G.bids[j][i]*(1-G.p[j]);
	            if(argmax<0 || ebid[argmax]<ebid[j])                 // first positive bid
	                argmax = j;
	        }
	    }
	    int j = argmax;             // i matches with j
	    M[i] = j;               // add the edge to matching
        if(j!=-1){
    	    u[i] = ebid[j];
    	    G.r[j] += G.bids[j][i]*G.p[j];  // update revenue of j 
    	    
    	    W += min(G.L[j], G.bids[j][i]);   // add to budget 
    	    Wf += (G.bids[j][i]-G.L[j]>0)?G.bids[j][i]-G.L[j] : 0;   // add to fake budget 
    	    G.L[j] -= min(G.L[j], G.bids[j][i]);   // update limit of j 
    	    algo+=G.bids[j][i];
        }
	}
	
	cout<<"Algorithm Matching-----\n";
	for(int i=0; i<n; i++)
	    cout<<"("<<i<<", "<<M[i]<<") "; 
	cout<<"\nAlgo="<<algo<<" W="<<W<<" Wf="<<Wf<<endl;

    vector<int> optimalMatch(n,-1);
    double OPT = opt(G, optimalMatch, 0);
    
	cout<<"Optimal Matching-----\n";
    for(int i=0; i<n; i++)
	    cout<<"("<<i<<", "<<optimalMatch[i]<<") "; 
    cout<<"\nOPT = "<<OPT<<endl;
    cout<<"Competitive Ratio: "<<W/OPT<<" vs (1-1/e)="<<1-exp(-1)<<endl<<"-----------------------"<<endl;
	return;
}

int main(){
    int t = 50;
    while(t--)
        test();
}