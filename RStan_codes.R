#Install and load the libraries
library("rstan")
library("rstudioapi")
#Base model
stanmodelcode_TN_base = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] av_prec;      //normalized precipitation
vector [nr] av_prec2;     //scaled precipitation
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag;                   //agriculture area in subwatershed
vector [nl] dev;               //Urban area in subwatershed
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs

}

transformed data{

}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d;        //Developed EC 
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 2> sigma_res;		//Model residual SD
real <lower=0, upper = 5> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads

}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] A;			//To compile Agriculture with PIC
vector [nr] D;			//To compile urban with PIC
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A;			//Agriculture vector
vector [nr] D_lc;		//Urban  vector		
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
vector [nr] ly_hat;
vector [nr] y;
// Loop to determine export for each watershed-year
for(i in 1:nr){
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc[i]= sum((ag[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc[i]= sum((dev[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

    //Adding precipitation impact coefficient to land export 

A[i]=  Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc[i];
D[i] = Be_d * pow(av_prec2[i],pic_p[2]) .* D_lc[i];
W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
  //Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(900,700);  //Prior for agriculture
Be_d ~ normal(800,300);  //Prior for development
Be_w ~ normal(200,200);   //Prior for undeveloped

Be_ch ~ normal(0.01,0.005);  //Prior for chickens
Be_h ~ normal(0.04,0.02);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.1);   //Prior for point source delivery
sigma_res ~ normal(0,20); //st error of the model
sigma_w ~ normal(0,300);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow


Sn ~ normal(.14,.05);    //Prior for stream retention rate
Sn2 ~ normal(11,2);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
ly_hat=log((y_hat /10000) + 10);   // vector to get log of TN load (y_hat)

ly ~ normal(ly_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
y=exp(ly)-10; 
load ~ normal(y,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'
#Buffer model
stanmodelcode_TN_buffer = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] av_prec;      //normalized precipitation
vector [nr] av_prec2;     //scaled precipitation
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag_buf;               //buffered agriculture area in subwatershed
vector [nl] ag_unbuf;             //unbuffered agriculture area in subwatershed
vector [nl] dev_buf;               //buffered urban area in subwatershed
vector [nl] dev_unbuf;            //Unbuffered urban area in subwatershed
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs

}

transformed data{

}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d;        //Developed EC 
real<lower =0, upper=1> mult;  //buffer efficiency (fraction) 
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 2> sigma_res;		//Model residual SD
real <lower=0, upper = 5> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads

}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] A;			//To compile Agriculture with PIC
vector [nr] D;			//To compile urban with PIC
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A_lc_buf;   //buffered agriculture vector
vector [nr] A_lc_unbuf; //Unbuffered agriculture vector
vector [nr] D_lc_buf;   //buffered urban vector		
vector [nr] D_lc_unbuf;   //Unbuffered urban vector		
vector [nr] D_buf;			//buffered urban vector
vector [nr] D_unbuf;			//buffered urban vector
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
vector [nr] ly_hat;
vector [nr] y;
// Loop to determine export for each watershed-year
for(i in 1:nr){
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc_buf[i]= sum((ag_buf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
A_lc_unbuf[i]= sum((ag_unbuf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_buf[i]= sum((dev_buf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_unbuf[i]= sum((dev_unbuf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

    //Adding precipitation impact coefficient to land export 

A[i]=  (1-mult)*Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc_buf[i]+Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc_unbuf[i];
D_buf[i] = (1-mult)*Be_d * pow(av_prec2[i],pic_p[2]) .* D_lc_buf[i];
D_unbuf[i] = Be_d* pow(av_prec2[i],pic_p[2]) .* D_lc_unbuf[i];
D[i]   = D_buf[i]+D_unbuf[i];
W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
  //Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(900,700);  //Prior for agriculture
Be_d ~ normal(800,300);  //Prior for development
mult~normal(0.45,0.25);   //Prior for buffer efficacy
Be_w ~ normal(200,200);   //Prior for undeveloped

Be_ch ~ normal(0.01,0.005);  //Prior for chickens
Be_h ~ normal(0.04,0.02);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.1);   //Prior for point source delivery
sigma_res ~ normal(0,20); //st error of the model
sigma_w ~ normal(0,300);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow


Sn ~ normal(.14,.05);    //Prior for stream retention rate
Sn2 ~ normal(11,2);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
ly_hat=log((y_hat /10000) + 10);   // vector to get log of TN load (y_hat)

ly ~ normal(ly_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
y=exp(ly)-10; 
load ~ normal(y,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'
#SCM model
stanmodelcode_TN_scm = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] av_prec;      //normalized precipitation
vector [nr] av_prec2;     //scaled precipitation
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag;               //agriculture area in subwatershed
vector [nl] dev_scm;               //urban area in subwatershed with SCM
vector [nl] dev_noscm;            //urban area in subwatershed without SCM
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs

}

transformed data{

}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d;        //Developed EC 
real<lower =0, upper=1> mult;  //scm efficiency (fraction) 
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 2> sigma_res;		//Model residual SD
real <lower=0, upper = 5> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads

}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] A;			//To compile Agriculture with PIC
vector [nr] D;			//To compile urban with PIC
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A_lc;   //buffered agriculture vector
vector [nr] D_lc_scm;   //urban (with SCM) vector		
vector [nr] D_lc_noscm;   //urban (without SCM) vector		
vector [nr] D_scm;			//urban (with SCM) vector
vector [nr] D_noscm;			//urban (without SCM) vector
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
vector [nr] ly_hat;
vector [nr] y;
// Loop to determine export for each watershed-year
for(i in 1:nr){
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc[i]= sum((ag[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_scm[i]= sum((dev_scm[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_noscm[i]= sum((dev_noscm[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

    //Adding precipitation impact coefficient to land export 

A[i]=  Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc[i];
D_scm[i] = (1-mult)*Be_d * pow(av_prec2[i],pic_p[2]) .* D_lc_scm[i];
D_noscm[i] = Be_d* pow(av_prec2[i],pic_p[2]) .* D_lc_noscm[i];
D[i]   = D_scm[i]+D_noscm[i];
W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
  //Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(900,700);  //Prior for agriculture
Be_d ~ normal(800,300);  //Prior for development
mult~normal(0.3,0.2);  //Prior for scm efficacy
Be_w ~ normal(200,200);   //Prior for undeveloped

Be_ch ~ normal(0.01,0.005);  //Prior for chickens
Be_h ~ normal(0.04,0.02);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.1);   //Prior for point source delivery
sigma_res ~ normal(0,20); //st error of the model
sigma_w ~ normal(0,300);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow


Sn ~ normal(.14,.05);    //Prior for stream retention rate
Sn2 ~ normal(11,2);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
ly_hat=log((y_hat /10000) + 10);   // vector to get log of TN load (y_hat)

ly ~ normal(ly_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
y=exp(ly)-10; 
load ~ normal(y,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'
#Buffer+SCM model
stanmodelcode_TN_scm_buffer = '

data { 
int nl;         //subwatersheds separated by year                            
int nr;         //incremental watershed-years
int nd;         // dischargers  
vector [nr] load;       //WRTDS load at LMS
vector [nr] SD;         	//WRTDS load SD
vector [nl] chick;        //number of chickens in subwatershed
vector [nl] cow;          //number of cows in subwatershed
vector [nl] hog;          //number of hogs (swine) in subwatershed
int wsd [nr];             //count variable for LMSs
int wshed_size;           //number of LMS watersheds
vector [nr] increm_area;  //Incremental area for each loading station
vector [nr] av_prec;      //normalized precipitation
vector [nr] av_prec2;     //scaled precipitation
vector [nr] up_t_load1;   //upstream loading for nested wsds
vector [nr] up_t_load2;   //upstream loading for nested wsds
vector [nr] up_t_load3;   //upstream loading for nested wsds
vector [nr] str_loss_load1;       //stream losses of upstream loading
vector [nr] str_loss_load2;       //stream losses of upstream loading
vector [nr] str_loss_load3;       //stream losses of upstream loading
vector [nr] res_loss_load1;       //reservoir losses of upstream loading
vector [nr] res_loss_load2;       //reservoir losses of upstream loading
vector [nr] res_loss_load3;       //reservoir losses of upstream loading
vector [nd] d_loss_str;           //stream losses of dischargers to LMSs
vector [nd] d_loss_res;           //reservoir losses of dischargers to LMSs
vector [nd] d_vals;               //discharger loadings
vector [nl] l_loss_str;           //stream losses for subwatersheds
vector [nl] l_loss_res;           //reservoir losses for subwatersheds
vector [nl] ag_buf;               //buffered agriculture area in subwatershed
vector [nl] ag_unbuf;             //unbuffered agriculture area in subwatershed
vector [nl] scm_buf;              //buffered urban area with SCM in subwatershed
vector [nl] scm_unbuf;            //unbuffered urban area with SCM in subwatershed
vector [nl] buf;                  //buffered urban area without SCM in subwatershed
vector [nl] unbuf;                //Unbuffered urban area without SCM in subwatershed
vector [nl] wild;			            //undeveloped area in subwatershed
vector [nl] tot_l;			          //total size of subwatershed
int l_start [nr];			            //index to link subwatersheds to LMSs
int l_end [nr];			              //index to link subwatersheds to LMSs
int d_start [nr];			            //index to link dischargers to LMSs
int d_end [nr];			              //index to link dischargers to LMSs

}

transformed data{

}

parameters {  

real<lower =0> Be_a;        //Agriculture EC
real<lower =0> Be_d;        //Developed EC 
real<lower =0, upper=1> mult1;  //buffer efficiency (fraction) 
real<lower =0, upper=1> mult2;  //scm efficiency (fraction) 
real<lower =0> Be_w;       //Undeveloped EC
real<lower =0, upper = 1> Be_ch;       //Chickens EC
real<lower =0, upper = 1> Be_h;       //Hogs (swine) EC
real<lower =0, upper = 1> Be_cw;       //Cows EC
real <lower =0, upper = 1> Sn;          // Stream retention
real <lower =1, upper = 60> Sn2;        //Reservoir retention
real <lower =0, upper = 0.40> PIC_q;      //PIC for retention
vector<lower =0, upper = 10> [7] pic_p;   //PIC for land classes
real <lower=0> Be_dch;      			//Point source DC 
real<lower=0, upper = 2> sigma_res;		//Model residual SD
real <lower=0, upper = 5> sigma_w;		//Random effect SD
vector [wshed_size] alpha;			//# of watersheds
real<lower = 0, upper = 2> sigma_B1;		//PIC SD
real<lower = 0, upper = 3> Bp_mean;		//PIC mean
vector [nr] ly;					// unknown true loads

}

transformed parameters {


}

model {
vector [nr] tot; 			// Sum of all loadings from all sources
vector [nr] sigma;		//SD for watershed random effects
vector [nr] y_hat;		//total loadings plus random effects minus waterbody losses
vector [nr] A;			//To compile Agriculture with PIC
vector [nr] D;			//To compile urban with PIC
vector [nr] W;			//To compile undeveloped urban with PIC
vector [nr] Dch; 		 	//To compile dischargers with stream/reservoir losses
vector [nr] alpha_vals;		// Watershed indicator
vector [nr] A_lc_buf;   //buffered agriculture vector
vector [nr] A_lc_unbuf; //Unbuffered agriculture vector
vector [nr] scm_lc_buf;   //buffered urban (with SCM) vector	
vector [nr] scm_lc_unbuf;   //unbuffered urban (with SCM) vector	
vector [nr] D_lc_buf;   //buffered urban (without SCM) vector		
vector [nr] D_lc_unbuf;   //Unbuffered urban (without SCM) vector		
vector [nr] D_buf;			//buffered urban vector
vector [nr] D_unbuf;			//unbuffered urban vector
vector [nr] W_lc;		//Undeveloped  vector
vector [nr] Disch;		//point source dischargers
vector [nr] C;   			//chickens for adding PIC
vector [nr] H;			//hogs for adding PIC
vector [nr] Cw;			//cows for adding PIC
vector [nr] C_r;   		//chickens for aggregating subwatersheds
vector [nr] H_r;			//swine for aggregating subwatersheds
vector [nr] Cw_r;		//cows  for aggregating subwatersheds
vector [nr] y;				//loading
vector [nr] ly_hat;
vector [nr] y;
// Loop to determine export for each watershed-year
for(i in 1:nr){
    //Looping to aggregate subwatersheds and corresponding stream and reservoir retention.
A_lc_buf[i]= sum((ag_buf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
A_lc_unbuf[i]= sum((ag_unbuf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

scm_lc_buf[i]= sum((scm_buf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
scm_lc_unbuf[i]= sum((scm_unbuf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

D_lc_buf[i]= sum((buf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
D_lc_unbuf[i]= sum((unbuf[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
W_lc[i]= sum((wild[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Disch[i]=sum((d_vals[d_start[i]:d_end[i]]) .* exp((-Sn/(1+PIC_q*av_prec[i])) * d_loss_str[d_start[i]:d_end[i]]) .* exp((-Sn2/(1+PIC_q*av_prec[i])) ./ d_loss_res[d_start[i]:d_end[i]]));
C_r[i]= sum((chick[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
H_r[i]= sum((hog[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));
Cw_r[i]= sum((cow[l_start[i]:l_end[i]]) .* (exp((-Sn/(1+PIC_q*av_prec[i])) * l_loss_str[l_start[i]:l_end[i]])) .* (exp((-Sn2/(1+PIC_q*av_prec[i])) ./l_loss_res[l_start[i]:l_end[i]])));

    //Adding precipitation impact coefficient to land export 

A[i]=  (1-mult1)*Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc_buf[i]+Be_a * pow(av_prec2[i],pic_p[1]) .* A_lc_unbuf[i];
D_buf[i] = (1-mult1)*(1-mult2)*Be_d * pow(av_prec2[i],pic_p[2]) .* scm_lc_buf[i]+(1-mult1) +Be_d * pow(av_prec2[i],pic_p[2]) .* D_lc_buf[i] ;
D_unbuf[i] = Be_d * pow(av_prec2[i],pic_p[2]) .* D_lc_unbuf[i]+(1-mult2)*Be_d * pow(av_prec2[i],pic_p[2]) .* scm_lc_unbuf[i];

D[i]   = D_buf[i]+D_unbuf[i];
W[i] = Be_w * pow(av_prec2[i],pic_p[4]) .* W_lc[i];

C[i] = Be_ch * pow(av_prec2[i],pic_p[5]) .* C_r[i];
H[i] = Be_h * pow(av_prec2[i],pic_p[6]) .* H_r[i];
Cw[i] = Be_cw * pow(av_prec2[i],pic_p[7]) .* Cw_r[i];

Dch[i] = Be_dch * Disch[i];
}


for (i in 1:nr){
  //Loop to determine random effect for each watershed
w= wsd[i];
sigma[i] = sqrt(pow(SD[i],2)+pow(sigma_res,2));   //with regular sigma values
alpha_vals[i] = alpha[w];

}
//Sum loadings from all sources
tot =  A + D + W + C + H + Cw + Dch;
//Add random effects to source loadings and subtract losses from upstream loads
y_hat = tot + alpha_vals - up_t_load1 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load1) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load1))) - up_t_load2 .* (1-(exp((-Sn ./ (1+PIC_q*av_prec)) .* str_loss_load2) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load2))) - up_t_load3 .* (1-(exp((-Sn ./(1+PIC_q*av_prec)) .* str_loss_load3) .* exp((-Sn2 ./ (1+PIC_q*av_prec)) ./ res_loss_load3)));


//priors
Be_a ~ normal(900,700);  //Prior for agriculture
Be_d ~ normal(800,300);  //Prior for development
mult1~normal(0.45,0.25);   //Prior for buffer efficacy
mult2~normal(0.3,0.2);   //Prior for buffer efficacy
Be_w ~ normal(200,200);   //Prior for undeveloped

Be_ch ~ normal(0.01,0.005);  //Prior for chickens
Be_h ~ normal(0.04,0.02);  //Prior for hogs (swine)
Be_cw ~ normal(0,5);  //Uninformed Prior for cows

Be_dch ~ normal(1,.1);   //Prior for point source delivery
sigma_res ~ normal(0,20); //st error of the model
sigma_w ~ normal(0,300);     //st. deviation of random effect hyperdistribution
alpha ~ normal(0,sigma_w);    //watershed random effects
sigma_B1 ~ normal(0,.5);  //st. deviation of precipitation coefficient
Bp_mean ~ normal(1,.5);    //mean PIC for hyperdistribution
pic_p[1] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for ag
pic_p[2] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for pre-1980 deve
pic_p[3] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for post-1980 dev
pic_p[4] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for undev
pic_p[5] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for chicken
pic_p[6] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for swine
pic_p[7] ~ normal(Bp_mean,sigma_B1);  //precipitation distribution for cow


Sn ~ normal(.14,.05);    //Prior for stream retention rate
Sn2 ~ normal(11,2);   //Prior for reservoir retention rate
PIC_q ~ normal(0,1);    //prior for PIC for retention
ly_hat=log((y_hat /10000) + 10);   // vector to get log of TN load (y_hat)

ly ~ normal(ly_hat,sigma_res);        //parameter that calibrates ly_hat with ly () 
y=exp(ly)-10; 
load ~ normal(y,SD);                // load = WRTDS estimate, SD = WRTDS sd

}

generated quantities {

}
'

data_base <- readRDS("./data_base.rds") #load input data set for the base model, if necessary
#run the base model (data is the list of data sets. For other parameters, return to the function description)
model_base = stan(model_code=stanmodelcode_TN_base, data=data_base, iter=iter, 
             warmup=warmup, thin=, chains=3,cores=3,
             control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))

data_buffer <- readRDS("./data_buffer") #load input data set for the base model, if necessary
#run the buffer model 
model_buffer = stan(model_code=stanmodelcode_TN_buffer, data=data_buffer, iter=iter, 
                  warmup=warmup, thin=, chains=3,cores=3,
                  control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))

data_scm <- readRDS("./data_scm") #load input data set for the base model, if necessary
#run the scm model 
model_scm = stan(model_code=stanmodelcode_TN_scm, data=data_scm, iter=iter, 
                    warmup=warmup, thin=, chains=3,cores=3,
                    control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))

data_buffer_scm <- readRDS("./data_buffer_scm") #load input data set for the base model, if necessary
#run the scm and buffer model 
model_buffer_scm = stan(model_code=stanmodelcode_TN_scm_buffer, data=data_buffer_scm, iter=iter, 
                 warmup=warmup, thin=, chains=3,cores=3,
                 control = list(adapt_delta =adapt_delta ,max_treedepth =max_treedepth ))