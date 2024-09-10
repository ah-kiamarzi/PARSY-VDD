PI_L1 float n;
PI_L1 int i2;
PI_L1 static int idx[NUM_CORES];

PI_L1 float zero = 0;
PI_L1 float one = 1;

PI_L1 static float buffer[NUM_CORES];
PI_L1 float rk = 0.0;
PI_L1 static float temp[NUM_CORES];

void GS_QR(float **Q, float **R, float **QRinput,int N_ROW,int N_COL);

#pragma GCC push_options
#pragma GCC optimize("-O3")
__attribute__((noinline)) void norm(float **R, float *v, int row, int column, int k, int blockSize_ROW, int start_ROW)
{
	int i2 = 0;
	float n = 0.0f;
	int j, idx;

	int core_id = pi_core_id();

	int start_matrice = start_ROW;

	buffer[core_id] = 0;
	idx = 0;

	for (j = start_ROW; (j + 1) < start_ROW + (blockSize_ROW); j += 2){
		float t0 = v[idx + start_matrice];
		float t1 = v[idx + 1 + start_matrice];
		buffer[core_id] = buffer[core_id] + t0 * t0;
		buffer[core_id] = buffer[core_id] + t1 * t1;
		idx += 2;
	}

	if (j < start_ROW + blockSize_ROW){
		buffer[core_id] = buffer[core_id] + v[idx + start_matrice] * v[idx + start_matrice];
	}
	
    #if NUM_CORES > 1
	pi_cl_team_barrier();
	BarrierCounter
    #endif

	if (core_id == 0){
		for (j = 0; j < NUM_CORES; j++){
			n += buffer[j];
		}
		R[k][k] = Sqrt(n);
		rk = 1 / R[k][k];
	}
    
    #if NUM_CORES > 1
	pi_cl_team_barrier();
	BarrierCounter
    #endif
}

void GS_QR(float **Q, float **R, float **QRinput, int N_ROW, int N_COL){
	int j;
	float RTEMP = 0;

	int core_id = pi_core_id();

	int blockSize_ROW = N_ROW / NUM_CORES;
	int start_ROW = core_id * blockSize_ROW;
	if (core_id == (NUM_CORES - 1))
	{
		blockSize_ROW = N_ROW - (NUM_CORES - 1) * blockSize_ROW;
	}
	int end_ROW = start_ROW + blockSize_ROW;

	for (int k = 0; k < N_COL; k++){
		for (j = start_ROW; ((j + 1) < end_ROW); j += 2){
			float in0 = QRinput[k][j];
			float in1 = QRinput[k][j + 1]; 
			Q[k][j] = in0;
			Q[k][j + 1] = in1;
		}
		if (j < end_ROW){
			float in0 = QRinput[k][j];
			Q[k][j] = in0;
		}

		#if NUM_CORES > 1
				pi_cl_team_barrier();
				BarrierCounter
		#endif

		for (int i = 0; i < k; i++){
			temp[core_id] = 0;
			for (j = start_ROW; ((j + 1) < end_ROW); j += 2){
				float ji0 = Q[i][j];
				float jk0 = Q[k][j];
				float ji1 = Q[i][j + 1];
				float jk1 = Q[k][j + 1];
				temp[core_id] = temp[core_id] + (ji0 * jk0) + (ji1 * jk1);
			}
			if ((j < end_ROW)){
				float ji0 = Q[i][j];
				float jk0 = Q[k][j];
				temp[core_id] = temp[core_id] + (ji0 * jk0);
			}
			
            #if NUM_CORES > 1
			pi_cl_team_barrier();
			BarrierCounter
            #endif

			#if NUM_CORES > 1
			RTEMP = temp[0] + temp[1];
			#else
			RTEMP = temp[0];
            #endif

            for (j = 2; j < NUM_CORES; j++){
				RTEMP += temp[j];
			}

			R[i][k] = RTEMP;

			for (j = start_ROW; (j + 1 < end_ROW); j += 2){
				float Rik = RTEMP;
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				float Qji1 = Q[i][j + 1];
				float Qjk1 = Q[k][j + 1];

				Q[k][j] = Qjk0 - Rik * Qji0;
				Q[k][j + 1] = Qjk1 - Rik * Qji1;
			}
			if (j < end_ROW){
				float Rik = RTEMP;
				float Qji0 = Q[i][j];
				float Qjk0 = Q[k][j];
				Q[k][j] = Qjk0 - Rik * Qji0;
			}
		}

		norm(R, &Q[k][0], N_ROW, N_ROW, k, blockSize_ROW, start_ROW);

		for (j = start_ROW; (j + 1) < (int)(start_ROW + (blockSize_ROW & 0xFFFFFFFE)); j += 2){
			float Qjk0 = Q[k][j];
			float Qjk1 = Q[k][j + 1];
			Q[k][j] = Qjk0 * rk;
			Q[k][j + 1] = Qjk1 * rk;
		}
		if (j < (int)(start_ROW + (blockSize_ROW & 0xFFFFFFFE))){
			float Qjk0 = Q[k][j];
			Q[k][j] = Qjk0 * rk;
		}

		if (blockSize_ROW & 0x0001){
			Q[k][j] = Q[k][j] * rk;
		}
	}
	return;
}
