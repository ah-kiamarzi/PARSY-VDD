#define MODEL_ORDER 16
#define NSMP 480

#define N 2048
#define NCHUNK 1 
#define NO (NSMP-MODEL_ORDER)
#define NO2 (NO-3*MODEL_ORDER/2)


#define MAX_ROWS 720


#if PPERF
	int numBarr = 0;
#endif

#include "pmsis.h"
#include "util.h"
#include "stdio.h"
#include "math.h"
#include "input.inc"
#include "stdlib.h"
#include "stdbool.h"
#include "cos2048_4.inc"

#if GS
	#include "GS_NEW.h"
#else
	#include "HH.h"
#endif

#if PPERF
	definePrefArr
#endif


PI_L1 float theta[MODEL_ORDER];
PI_L1 float PSD_RES[N];
PI_L1 float sigma[NUM_CORES];
PI_L1 float Yi[MODEL_ORDER];


PI_L1 float u[NO+MODEL_ORDER];

float **Q, **R, **phi_transposed;

// EXTRA VARIABLE FOR HH
// float **Q_temp, **R_temp, **v_temp;



void cluster_main();

void pe_entry(void *arg){
	cluster_main();
}

void cluster_entry(void *arg){
	pi_cl_team_fork((NUM_CORES), pe_entry, 0);
}

static int test_entry(){

	// pi_perf_conf(1 << PI_PERF_CYCLES);

	struct pi_device cluster_dev;
	struct pi_cluster_conf cl_conf;
	struct pi_cluster_task cl_task;

	pi_cluster_conf_init(&cl_conf);

	cl_conf.id = 0;
	// cl_conf.icache_conf = PI_CLUSTER_MASTER_CORE_ICACHE_ENABLE |
	// 					// Enable the prefetch for all the cores, it's a 9bits mask (from bit 2 to bit 10), each bit correspond to 1 core
	// 					PI_CLUSTER_ICACHE_PREFETCH_ENABLE |
	// 					// Enable the icache for all the cores
	// 					PI_CLUSTER_ICACHE_ENABLE;
	pi_open_from_conf(&cluster_dev, &cl_conf);
	if (pi_cluster_open(&cluster_dev)){
		return -1;
	}

	pi_cluster_task(&cl_task, cluster_entry, NULL);

	pi_cluster_send_task_to_cl(&cluster_dev, &cl_task);

	pi_cluster_close(&cluster_dev);

	return 0;
}


int main()
{
// Set voltage to 0.8 and frequency to 370 MHz
	pi_pmu_voltage_set(PI_PMU_DOMAIN_CL, 800);
	pi_freq_set(PI_FREQ_DOMAIN_FC, 370*1000*1000);
	pi_freq_set(PI_FREQ_DOMAIN_CL, 370*1000*1000);
	pi_freq_set(PI_FREQ_DOMAIN_PERIPH, 370*1000*1000);


	return pmsis_kickoff((void *)test_kickoff);

    int ret = test_entry();
    pmsis_exit(ret);
	return 0;

}
void SI();


void memoryAllocation(){
	Q = (float **) pi_l1_malloc(0,MODEL_ORDER * sizeof(float *));
	for(int i = 0; i < MODEL_ORDER; i++){
		Q[i] = (float *) pi_l1_malloc(0,NO * sizeof(float));
	}
	phi_transposed = (float **) pi_l1_malloc(0,MODEL_ORDER * sizeof(float *));
	for(int i = 0; i < MODEL_ORDER; i++){
		phi_transposed[i] = (float *) pi_l1_malloc(0,NO * sizeof(float));
	}
	#if GS
		R = (float **) pi_l1_malloc(0,MODEL_ORDER * sizeof(float *));
		for(int i = 0; i < MODEL_ORDER; i++){
			R[i] = (float *) pi_l1_malloc(0,MODEL_ORDER * sizeof(float));
		}
	#else
		R = (float **) pi_l2_malloc(MODEL_ORDER * sizeof(float *));
		for(int i = 0; i < MODEL_ORDER; i++){
			R[i] = (float *) pi_l2_malloc(NO * sizeof(float));
		}
		Q_temp = (float **) pi_l2_malloc(MODEL_ORDER * sizeof(float *));
		for(int i = 0; i < MODEL_ORDER; i++){
			Q_temp[i] = (float *) pi_l2_malloc(NO * sizeof(float));
		}
		R_temp = (float **) pi_l2_malloc(MODEL_ORDER * sizeof(float *));
		for(int i = 0; i < MODEL_ORDER; i++){
			R_temp[i] = (float *) pi_l2_malloc(NO * sizeof(float));
		}

		v_temp = (float **) pi_l2_malloc(NO * sizeof(float *));
		for(int i = 0; i < NO; i++){
			v_temp[i] = (float *) pi_l2_malloc(NO * sizeof(float));
		}

	#endif
}
void dataPreparation(int thisFNO, bool mode, int core_id){//parallelized

	float input_load_temp[MODEL_ORDER];
	#if NUM_CORES > 1
	int blockSize = thisFNO/NUM_CORES;
	int start = core_id*blockSize;

	if(core_id==(NUM_CORES - 1)){
		blockSize = thisFNO - (NUM_CORES - 1)* blockSize;}
	int end = blockSize + start; 

	#endif

	if(mode){ 
		#if NUM_CORES > 1
		for(int i=(start + MODEL_ORDER/2);i<(end + MODEL_ORDER/2);i++){
			for(int j=0;j<MODEL_ORDER/2;j++){
				phi_transposed[j][i-MODEL_ORDER/2] = -input[i-j+MODEL_ORDER-1];
				phi_transposed[j+MODEL_ORDER/2][i-MODEL_ORDER/2] =  u[i-j+MODEL_ORDER-1];					
			}
		}
		pi_cl_team_barrier();
		BarrierCounter

		#else
		for(int i=MODEL_ORDER/2;i<NO2+MODEL_ORDER/2;i++){
			for(int j=0;j<MODEL_ORDER/2;j++){
				phi_transposed[j][i-MODEL_ORDER/2] = -input[i-j+MODEL_ORDER-1];
				phi_transposed[j+MODEL_ORDER/2][i-MODEL_ORDER/2] =  u[i-j+MODEL_ORDER-1];
			}
		}
		#endif
	}
	#if NUM_CORES > 1
		pi_cl_team_barrier();
		BarrierCounter
	#endif
	

	#if NUM_CORES > 1

	for(int j = 0; j < MODEL_ORDER; j++){
		input_load_temp[j] = -input[start+j];
	}	

	for(int i = start; i < end; i+=1){
		for(int j = 0; j < MODEL_ORDER; j++){
			if(!mode){
					phi_transposed[j][i] = input_load_temp[((MODEL_ORDER-1-j)+(i-start))%MODEL_ORDER];
			}
		}
		if(!mode){
			input_load_temp[(i-start)%MODEL_ORDER] = -input[i+MODEL_ORDER];
		}
	}
	pi_cl_team_barrier();
	BarrierCounter
	#else
	
	for(int i = 0; i < thisFNO; i+=1){
		for(int j = 0; j < MODEL_ORDER; j++){
			if(!mode){
				phi_transposed[j][i] = input_load_temp[((MODEL_ORDER-1-j)+i)%MODEL_ORDER];
			}
		}
		if(!mode){
			input_load_temp[i%MODEL_ORDER] = -input[i+MODEL_ORDER];
		}
	}	
	#endif

}
void VCUpdate(int thisFNO, bool mode,int core_id){//vectore coeffecient update//parallelized
	#if NUM_CORES > 1
	int blockSize = MODEL_ORDER/NUM_CORES;
	int start = core_id*blockSize;

	if(core_id==(NUM_CORES - 1)){
		blockSize = MODEL_ORDER - (NUM_CORES - 1)* blockSize;}
	int end = blockSize + start;

	for(int i=start; i < end; i++){
		Yi[i] = 0;
		for(int w=0; w<thisFNO; w++){
			if(mode){
				Yi[i] += Q[i][w] * input[w+3*MODEL_ORDER/2];
			}else{
				Yi[i] += Q[i][w] * input[w+MODEL_ORDER];
			}
		}
	}
	pi_cl_team_barrier();
	BarrierCounter
	#else
	for(int i=0; i < MODEL_ORDER; i++){
		Yi[i] = 0;
		for(int w=0; w<thisFNO; w++){
			if(mode){
				Yi[i] += Q[i][w] * input[w+3*MODEL_ORDER/2];
			}else{
				Yi[i] += Q[i][w] * input[w+MODEL_ORDER];
			}
		}
	}	
	#endif
}

void sysParEstimation(){//SysId parameter estimation//sequential
	float sum;
	for(int k=MODEL_ORDER-1; k>-1; k--){
		if(k==MODEL_ORDER-1)
			theta[k] = Yi[k]/R[k][k];
		else{
			sum = 0;
			for(int i=k+1; i<MODEL_ORDER;i++){
				#if GS
					sum+=theta[i]*R[k][i];				
				#else
					sum+=theta[i]*R[i][k];
				#endif
			}
			theta[k]=(Yi[k]-sum)/R[k][k];
		}
	}
}

void noiseVarEstimation(int thisFNO, bool mode, int core_id){//noise variance estimation//parallelized
	
	#if NUM_CORES > 1
	int blockSize = (thisFNO-MODEL_ORDER)/NUM_CORES;
	int start = core_id*blockSize+MODEL_ORDER;

	if(core_id==(NUM_CORES - 1)){
		blockSize = (thisFNO-MODEL_ORDER) - (NUM_CORES - 1)* blockSize;}
	int end = blockSize + start;

	sigma[core_id] = 0;
	for(int j=start; j<end; j++){
		float sum_theta_y = 0;
		for(int i=0; i<MODEL_ORDER; i++){
			if(mode){
				sum_theta_y += theta[i] * (-input[MODEL_ORDER + j-i-1]);
			}else{
				sum_theta_y += theta[i] * (-input[j-i-1]);
			}
		}
		if(mode){
			u[j] = input[j+MODEL_ORDER] - sum_theta_y;
		}else{
			u[j] = input[j] - sum_theta_y;
		}

		sigma[core_id] += u[j] * u[j];
	}
	sigma[core_id] /= (thisFNO-MODEL_ORDER-1);

	pi_cl_team_barrier();
	BarrierCounter
	if(core_id == 0){
		for (int i = 1; i < NUM_CORES; i++){
			sigma[0] += sigma[i];
		}
	}
	pi_cl_team_barrier();
	BarrierCounter
	#else
	sigma[core_id] = 0;
	for(int j=MODEL_ORDER; j<thisFNO; j++){
		float sum_theta_y = 0;
		for(int i=0; i<MODEL_ORDER; i++){
			if(mode){
				sum_theta_y += theta[i] * (-input[MODEL_ORDER + j-i-1]);
			}else{
				sum_theta_y += theta[i] * (-input[j-i-1]);
			}
		}
		if(mode){
			u[j] = input[j+MODEL_ORDER] - sum_theta_y;
		}else{
			u[j] = input[j] - sum_theta_y;
		}

		sigma[core_id] += u[j] * u[j];
	}
	sigma[core_id] /= (thisFNO-MODEL_ORDER-1);
	#endif
}

void QR_decompositon(int thisFNO,int core_id){
	#if GS
		initialMatrixConstPP(Q, MODEL_ORDER, thisFNO, zero, core_id);
		initialMatrixConstPP(R, MODEL_ORDER, MODEL_ORDER, zero, core_id);

		pi_perf_fc_reset();
		pi_perf_fc_start();

		GS_QR(Q, R, phi_transposed, thisFNO/NCHUNK, MODEL_ORDER);
		pi_perf_fc_stop();

		#if PPERF
			printPerfCounters
			if(pi_core_id() == 0){
				printf("NUMBARR = %d\n",numBarr);
			}
		#endif

		printMatrixPP(Q,thisFNO,MODEL_ORDER,"Q",1);

		if(core_id == 0){
			printMatrixPP(R,MODEL_ORDER,MODEL_ORDER,"R",0);
		}
		pi_cl_team_barrier();

	#else
		initialMatrixEYEPP(Q,MODEL_ORDER,thisFNO/NCHUNK);
		initialMatrixMatrixPP(R,MODEL_ORDER,thisFNO/NCHUNK,phi_transposed);
		initialMatrixConstPP(R_temp,MODEL_ORDER,thisFNO/NCHUNK,0,core_id);	
		HH_QR(Q, R, Q_temp, R_temp, v_temp, phi_transposed, thisFNO/NCHUNK, MODEL_ORDER);
	#endif

}

void AR_Func(int core_id){


		VCUpdate(NO,false,core_id);
		if(core_id == 0){
			sysParEstimation();
		}
		pi_cl_team_barrier();
		BarrierCounter

		noiseVarEstimation(NO,false,core_id);

		#if DEBUG
			if(core_id == 0){
				printf("AR SIGMA = %f\n",sigma[0]);
			}
			pi_cl_team_barrier();
		#endif
}

void ARMA_Func(int core_id){

	VCUpdate(NO2,true,core_id);
	if(core_id == 0){
		sysParEstimation();
	}
	pi_cl_team_barrier();
	BarrierCounter

	noiseVarEstimation(NO2,true,core_id);

	#if DEBUG
		if(core_id == 0){
			printf("ARMA SIGMA = %f\n",sigma[0]);
		}
		pi_cl_team_barrier();
	#endif
}


void PSD(int core_id){
	#if NUM_CORES > 1
	int blockSize = N/NUM_CORES;
	int start = core_id*blockSize;

	if(core_id==(NUM_CORES - 1)){
		blockSize = N - (NUM_CORES - 1)* blockSize;}
	int end = blockSize + start;
	#endif
	float A;
	float B;
	float C;
	float D;

	if(!ARMA){

		#if NUM_CORES > 1
		for (int ii = start; ii < end; ii++){
			A = one;
			B = zero;
			C = one;
			D = zero;

			for(int j = 0; j < MODEL_ORDER; j++){

				int index = (ii*j+ii) & 0x000007FF;
				float cosVal, sinVal;
				float thetaVal = theta[j];
				if(index < N/4){
					cosVal = cosLUT[index];
					sinVal = cosLUT[N/4 - index];
				}else if(index < N/2){
					cosVal = -cosLUT[N/2 - index];
					sinVal = cosLUT[index - N/4];
				}else if(index < 3*N/4){
					cosVal = -cosLUT[index - N/2];
					sinVal = -cosLUT[3*N/4 - index];
				}else{
					cosVal = cosLUT[N - index];
					sinVal = -cosLUT[index-3*N/4];
				}
				C += thetaVal*cosVal;
				D += thetaVal*sinVal;
			}
			float C2 = C*C;
			float D2 = D*D;
			PSD_RES[ii] =	(sigma[0]*sigma[0])/((C2 + D2)*(C2 + D2));
		}
		#else
		for (int ii = 0; ii < N; ii++){
			A = one;
			B = zero;
			C = one;
			D = zero;

			for(int j = 0; j < MODEL_ORDER; j++){

				int index = (ii*(j+1)) & 0x000007FF;
				float cosVal, sinVal;
				float thetaVal = theta[j];
				if(index < N/4){
					cosVal = cosLUT[index];
					sinVal = cosLUT[N/4 - index];
				}else if(index < N/2){
					cosVal = -cosLUT[N/2 - index];
					sinVal = cosLUT[index - N/4];
				}else if(index < 3*N/4){
					cosVal = -cosLUT[index - N/2];
					sinVal = -cosLUT[3*N/4 - index];
				}else{
					cosVal = cosLUT[N - index];
					sinVal = -cosLUT[index-3*N/4];
				}
				C += thetaVal*cosVal;
				D += thetaVal*sinVal;
			}
			float C2 = C*C;
			float D2 = D*D;
			PSD_RES[ii] =	(sigma[0]*sigma[0])/((C2 + D2)*(C2 + D2));
		}
		#endif
	}else{
		#if NUM_CORES > 1
		for (int ii = start; ii < end; ii++){
			A = one;
			B = zero;
			C = one;
			D = zero;
			for(int j = 0; j < MODEL_ORDER/2; j++){
				int index = (ii*j+ii) & 0x000007FF;
				float cosVal, sinVal;
				float thetaVal = theta[j];
				float thetaVal2 = theta[j+MODEL_ORDER/2];

				if(index < N/4){
					cosVal = cosLUT[index];
					sinVal = cosLUT[N/4 - index];
				}else if(index < N/2){
					cosVal = -cosLUT[N/2 - index];
					sinVal = cosLUT[index - N/4];
				}else if(index < 3*N/4){
					cosVal = -cosLUT[index - N/2];
					sinVal = -cosLUT[3*N/4 - index];
				}else{
					cosVal = cosLUT[N - index];
					sinVal = -cosLUT[index-3*N/4];
				}

				A += thetaVal*cosVal;
				B += thetaVal*sinVal;
				C += thetaVal2*cosVal;
				D += thetaVal2*sinVal;				

			}
			float C2 = C*C;
			float D2 = D*D;
			float temp1 = (A*C+B*D)*(A*C+B*D);
			float temp2 = (-A*D+B*C)*(-A*D+B*C);
			PSD_RES[ii] =	((sigma[0]*sigma[0])*((C2 + D2)*(C2 + D2)))/(temp1+temp2);
		}
		#else
		for (int ii = 0; ii < N; ii++){
			A = one;
			B = zero;
			C = one;
			D = zero;
			for(int j = 0; j < MODEL_ORDER/2; j++){
				int index = (ii*(j+1)) & 0x000007FF;
				float cosVal, sinVal;
				float thetaVal = theta[j];
				float thetaVal2 = theta[j+MODEL_ORDER/2];

				if(index < N/4){
					cosVal = cosLUT[index];
					sinVal = cosLUT[N/4 - index];
				}else if(index < N/2){
					cosVal = -cosLUT[N/2 - index];
					sinVal = cosLUT[index - N/4];
				}else if(index < 3*N/4){
					cosVal = -cosLUT[index - N/2];
					sinVal = -cosLUT[3*N/4 - index];
				}else{
					cosVal = cosLUT[N - index];
					sinVal = -cosLUT[index-3*N/4];
				}

				A += thetaVal*cosVal;
				B += thetaVal*sinVal;
				C += thetaVal2*cosVal;
				D += thetaVal2*sinVal;				

			}
			float C2 = C*C;
			float D2 = D*D;
			float temp1 = (A*C+B*D)*(A*C+B*D);
			float temp2 = (-A*D+B*C)*(-A*D+B*C);
			PSD_RES[ii] =	((sigma[0]*sigma[0])*((C2 + D2)*(C2 + D2)))/(temp1+temp2);
		}
		#endif
	}
}



void cluster_main()
{

	// pi_perf_conf(
	// (1<<PI_PERF_CYCLES) | 
	// (1<<PI_PERF_ACTIVE_CYCLES) | 
	// (1<<PI_PERF_INSTR)  | 
	// (1<<PI_PERF_LD_EXT) | 
	// (1<<PI_PERF_TCDM_CONT) | 
	// (1<<PI_PERF_LD_STALL) | 
	// (1<<PI_PERF_IMISS) | 
	// (1<<PI_PERF_LD_EXT_CYC) | 
	// (1<<PI_PERF_ST_EXT_CYC) | 
	// (1<<0x11) | 
	// (1<<0x12) |
	// (1<<0x13) | 
	// (1<<0x14));
	
	// pi_perf_conf((1<<PI_PERF_CYCLES) | (1<<PI_PERF_INSTR));



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//									 ASSIGN THE INITIAL VALUES HERE
//////////////////////////////////////////////////////////////////////////////////////////////////////////
	pi_perf_reset();
	pi_perf_start();
		SI();
	pi_perf_stop();

	#if PPERF
		printPerfCounters
		if(pi_core_id() == 0){
			printf("NUMBARR = %d\n",numBarr);
		}
	#endif
	pi_cl_team_barrier();




}

void SI(){//SysID

	int core_id = pi_core_id();

/************************************************/
	if(core_id == 0){
		memoryAllocation();
	}
	#if NUM_CORES > 1
	pi_cl_team_barrier();
	BarrierCounter
	#endif
/************************************************/




	dataPreparation(NO,false,core_id);
	if(ARMA){
		dataPreparation(NO2,true,core_id);
		QR_decompositon(NO2,core_id);
		ARMA_Func(core_id);
	}
	#if NUM_CORES>1
		pi_cl_team_barrier();
		BarrierCounter
	#endif

}
