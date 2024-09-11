PI_L1 float zero = 0;
PI_L1 float one = 1;
PI_L1 float two = 2;
PI_L1 float b[MODEL_ORDER];


PI_L1 float v[NO];
// PI_L2 float v[NO];


void HH_QR(float **Q, float **R, float **Q_temp, float **R_temp,float **v_temp, float **input,int N_ROW, int N_COL);


void house_opt(float *v,float *b, float *x, int n){
	float sigma = 0;
	float x0 = x[0];
	if(n == 1){
		sigma = 0;
	}else{
		int i;
		for (i = 1; i+1 < n; i+=2){
			float x0 = x[i];
			float x1 = x[i+1];
			float x0p2 = x0*x0;
			float x1p2 = x1*x1;
			sigma = sigma + x0p2+x1p2;
		}
		if(i < n){
			float x0 = x[i];
			float x0p2 = x0*x0;
			sigma = sigma + x0p2;
		}
	}
	v[0] = one;
	int i;
	for (i = 1; i+1 < n; i+=2){
		float x0 = x[i];
		float x1 = x[i+1];
		v[i] = x0;
		v[i+1] = x1;
	}
	if(i < n){
		float x0 = x[i];
		v[i] = x0;
	}

	if(sigma == 0){
		*b = 0;
	}else{
		float mu = Sqrt(x0*x0 + sigma);
		if(x0 <= 0){
			v[0] = x0 - mu;
		}else{
			v[0] = -sigma/(x0+mu);
		}
	}
	float v0 = v[0]; 
	float v0p2 = v0*v0;
	*b = two*v0p2/(sigma+v0p2);
	v[0] = one;
	for (int i = 1; i < n; i++){
		v[i] = v[i]/v0;
	}
}

void HH_QR(float **Q, float **R, float **Q_temp, float **R_temp,float **v_temp, float **input,int N_ROW, int N_COL){
	int core_id = pi_core_id();
	int l;
	float bj;
	#if NUM_CORES > 1
		int blockSize_ROW = N_ROW/NUM_CORES;
		int start_ROW = core_id*blockSize_ROW;

		if(core_id==(NUM_CORES - 1)){
			blockSize_ROW = N_ROW - (NUM_CORES - 1)* blockSize_ROW;}
		int end_ROW = start_ROW+blockSize_ROW;

		int blockSize_COL = N_COL/NUM_CORES;
		int start_COL = core_id*blockSize_COL;

		if(core_id==(NUM_CORES - 1)){
			blockSize_COL = N_COL - (NUM_CORES - 1)* blockSize_COL;}
		int end_COL = start_ROW+blockSize_COL;
		int start = 0;

	#endif

	for (int j = 0; j < N_COL; j++){
		if(core_id == 0){
			house_opt(&v[0],&b[j], &R[j][j], N_ROW-j);
		}
		#if NUM_CORES > 1
			pi_cl_team_barrier();
			BarrierCounter
		#endif
		bj = b[j];
		#if NUM_CORES > 1
			for (int i = start_ROW; ((i < N_ROW-j) && (i < end_ROW)); i++){
				float vi = v[i];
				for (int k = 0; k < N_ROW-j; k++){
					if(i == k){
						v_temp[i][k] =	one - bj * vi*vi;
					}else{
						v_temp[i][k] = -bj * vi*v[k];
					}
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for (int i = 0; i < N_ROW-j; i++){
				float vi = v[i];
				for (int k = 0; k < N_ROW-j; k++){
					if(i == k){
						v_temp[i][k] =	one - bj * vi*vi;
					}else{
						v_temp[i][k] = -bj * vi*v[k];
					}
				}
			}
		#endif

		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}
			for(int i = start; i < end_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;

					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i-j][l-j]*R[k][l];
					}
					R_temp[k][i] = temp;
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i-j][l-j]*R[k][l];
					}
					R_temp[k][i] = temp;
				}
			}
		#endif

		#if NUM_CORES > 1
			if(start_COL >= j){
				start = start_COL;
			}else if(end_COL >= j){
				start = j;
			}else{
				start = end_COL+1;
			}

			for(int i = start; i < end_COL; i++){
				int k;
				for (k = j; k+1 < N_ROW; k+=2){
					float R_temp_ik0 = R_temp[i][k];
					float R_temp_ik1 = R_temp[i][k+1];
					R[i][k] = R_temp_ik0;
					R[i][k+1] = R_temp_ik1;
				}
				if(k < N_ROW){
					float R_temp_ik0 = R_temp[i][k];
					R[i][k] = R_temp_ik0;
				}
			}

		#else
			for(int i = j; i < N_COL; i++){
				for (int k = j; k < N_ROW; k++){
					R[i][k] = R_temp[i][k];
				}
			}
		#endif

		if(j < N_ROW){
		#if NUM_CORES > 1
			if(start_ROW >= j+1){
				start = start_ROW;
			}else if(end_ROW >= j+1){
				start = j+1;
			}else{
				start = end_ROW+1;
			}

			int i;
			for(i = start; i+1 < end_ROW; i+=2){
				float v0 = v[1+i-j-1]; 
				float v1 = v[1+i+1-j-1]; 
				R[j][i] = v0;
				R[j][i+1] = v1;
			}
			if(i < end_ROW){
				float v0 = v[1+i-j-1]; 
				R[j][i] = v0;
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j+1; i < N_ROW; i++){
				R[j][i] = v[1+i-j-1];
			}

		#endif
		}
	}


	for (int j = N_COL-1; j >= 0; j--){
		v[j] = one;
		#if NUM_CORES > 1
			if(start_ROW >= j+1){
				start = start_ROW;
			}else if(end_ROW >= j+1){
				start = j+1;
			}else{
				start = end_ROW+1;
			}

			for(int i = start; i < end_ROW; i++){
				v[i] = R[j][i];
			}
		#else
			for(int i = j+1; i < N_ROW; i++){
				v[i] = R[j][i];
			}
		#endif


		#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}

			float bj = b[j];
			for (int i = start; i < end_ROW; i++){
				float vi = v[i];
				for (int k = j; k < N_ROW; k++){
					if(i == k){
						v_temp[i][k] = one - bj * vi*vi;
					}else{
						v_temp[i][k] =   - bj * vi*v[k];
					}
				}
			}
		#else
			float bj = b[j];
			for (int i = j; i < N_ROW; i++){
				float vi = v[i];
				for (int k = j; k < N_ROW; k++){
					if(i == k){
						v_temp[i][k] = one - bj * vi*vi;
					}else{
						v_temp[i][k] =   - bj * vi*v[k];
					}
				}
			}
		#endif

				#if NUM_CORES > 1
			if(start_ROW >= j){
				start = start_ROW;
			}else if(end_ROW >= j){
				start = j;
			}else{
				start = end_ROW+1;
			}
			for(int i = start; i < end_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i][l]*Q[k][l];
					}
					Q_temp[k][i] = temp;
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_ROW; i++){
				for (int k = j; k < N_COL; k++){
					float temp = 0;
					for (int l = j; l < N_ROW; l++){
						temp += v_temp[i][l]*Q[k][l];
					}
					Q_temp[k][i] = temp;
				}
			}
		#endif
	
		#if NUM_CORES > 1
			if(start_COL >= j){
				start = start_COL;
			}else if(end_COL >= j){
				start = j;
			}else{
				start = end_COL+1;
			}			
			for(int i = start; i < end_COL; i++){
				for (int k = j; k < N_ROW; k++){
					Q[i][k] = Q_temp[i][k];
				}
			}
			pi_cl_team_barrier();
			BarrierCounter
		#else
			for(int i = j; i < N_COL; i++){
				for (int k = j; k < N_ROW; k++){
					Q[i][k] = Q_temp[i][k];
				}
			}		
		#endif
	}
	return;
}