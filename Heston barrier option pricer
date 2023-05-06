#include <iostream>
#include <random>
#include <cmath>
using namespace std;
   //payout function at terminal date
double payout(double ST,char callput,double strike){
        if(callput=='c'){
            return max(ST-strike,0.0);
        }
        if(callput=='p'){
            return max(-ST+strike,0.0);
        }

    }
int main()
{
    // starting normal number generator
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,1.0);
    string readme="inputs:\ncallput - 'c' for call option and 'p' for put\ninout - 'i' for in option 'o' for out option\nB - barrier\nK - strike\nS0 - underlying price at time 0\nT - terminal date\nr - domestic rate\nq - foregin rate\nv0 - inintal square root of volatility in martingale mesure\nVT - long time avrege square root of volatility in martingale mesure\nkappa - mean reversion rate of volatility in martingale mesure\nrho - corelation of wiener processes in martingale mesure\nksi - volatility of volatility in martingale mesure\nM- number of monte carlo simulations\nN- number of steps in each simulation\n\nThe stochastic pocess of price in martingale mesure is given as:\ndS_t=(r-q)S_t dt + sqrt(v_t) dW1_t\ndv_t=kappa(VT-v_t) dt + ksi sqrt(v_t) d W2_t\nwhere corr(W1,W2)=rho\ninintial parameters are: S_0=S0, v_0=v0\n";




    double sample=distribution(generator);
    double B,K,S0,T,r,q,v0,VT,kappa,rho,ksi;
    int M,N;
    char callput,inout;
    cout<<readme<<endl<<endl<<"Enter callput,inout,B,K,S0,T,r,q,v0,VT,kappa,rho,ksi,M,N"<<endl;
    cin>>callput>>inout>>B>>K>>S0>>T>>r>>q>>v0>>VT>>kappa>>rho>>ksi>>M>>N;
    /* input example
    c o 0 500 553.82 0.06746032 0.05 0 0.1147304 0.1678504 0.1913857 -0.01523716 1.66512 10000 10000
    true value is 58.40603
    */
    double dt=T/(double) N;
    double P[M]; //array to hold monte carlo payouts
    for(int m=0;m<M;m++){
        double St=S0; double vt=v0;
        bool IN=0; //flag to remember if the barrier was breached for in options
        for(int n=0;n<N;n++){
            double dW1=distribution(generator)*sqrt(dt);

            double dW2=rho*dW1+sqrt(1-pow(rho,2))*distribution(generator)*sqrt(dt);
            St+=(r-q)*St*dt+sqrt(abs(vt))*St*dW1;
            //up out break
            if(S0<B & inout=='o' & St>B){
                    P[m]=0;
                    break;
            }
            //down out break
            if(S0>=B & inout=='o' & St<B){
                    P[M]=0;
                    break;
            }
            // down in flag
            if(S0>=B & inout=='i' & St<B){
                    IN=1;
            }
            // up in flag
            if(S0<B & inout=='i' & St>B){
                    IN=1;
            }

            vt+=kappa*(VT-vt)*dt+ksi*sqrt(abs(vt))*dW2;
        }
        if(inout=='o'){
            P[m]=payout(St,callput,K);
        }
        if(inout=='i'){
            if(IN==1){
                P[m]=payout(St,callput,K);
            }
            else{
                P[m]=0;
            }
        }
    }
    double price=0.0;
    double price2=0.0;
    for(int m=0;m<M;m++){
        price+=P[m];
        price2+=pow(P[m],2);
    }
    price=exp(-r*T)*price/(double)M;
    cout<<"Monte Carlo Price:\t\t"<<(double)price<<endl;
    double variance=((price2/(double)M-pow(price,2))*exp(-2*r*T))/(double)M;
    double alpha=sqrt(variance)*1.64;
    cout<<"90% confidence interval:\t"<<price-alpha<<"   "<<price+alpha;
    return 0;
}
