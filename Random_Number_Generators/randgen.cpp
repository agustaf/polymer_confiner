
#define W 32
#define R 32
#define M1 3
#define M2 24
#define M3 10

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define Identity(v) (v)

#define V0            STATE[state_i                   ]
#define VM1           STATE[(state_i+M1) & 0x0000001fU]
#define VM2           STATE[(state_i+M2) & 0x0000001fU]
#define VM3           STATE[(state_i+M3) & 0x0000001fU]
#define VRm1          STATE[(state_i+31) & 0x0000001fU]
#define newV0         STATE[(state_i+31) & 0x0000001fU]
#define newV1         STATE[state_i                   ]

#define FACT 2.32830643653869628906e-10

static unsigned int state_i = 0;
static unsigned int STATE[R];
static unsigned int z0, z1, z2;

void randgenbase::init_randgen(unsigned int* init){
	int j;
	state_i = 0;
	for (j = 0; j < R; ++j){
		STATE[j] = init[j];
	}
	return;
}
double randgenbase::randgen(){
	z0    = VRm1;
	z1    = Identity(V0)       ^ MAT0POS (8, VM1);
	z2    = MAT0NEG (-19, VM2) ^ MAT0NEG(-14,VM3);
	newV1 = z1                 ^ z2; 
	newV0 = MAT0NEG (-11,z0)   ^ MAT0NEG(-7,z1)    ^ MAT0NEG(-13,z2) ;
	state_i = (state_i + 31) & 0x0000001fU;
	return ((double) STATE[state_i]  * FACT);
}
void randgenbase::autoinit_randgen(){
	const int seed = time(NULL);
	srand(seed);
	const unsigned int max_size = std::numeric_limits<unsigned int>::max();
	unsigned int* firstvalues = 0;
	firstvalues = new unsigned int[R];
	for(int i=0; i<R; ++i){
		firstvalues[i] = (unsigned int)(double(max_size)* \
		double(rand())/double(RAND_MAX));
	}
	init_randgen(firstvalues);
	delete[] firstvalues;
	firstvalues = 0;
	return;
}
void randgenbase::advance_randgen(const int count){
	for(int i=0; i<count; ++i){
		randgen();
	}
}
static double randgenbase::gaussrandgen(){
	static int value_is_cached = 0;
	static double cached_value = 0.0;
	if(value_is_cached == 1){
		value_is_cached = 0;
		return(cached_value);
	}
	else{
		double randnum1 = 2.0*randgen() - 1.0;
		double randnum2 = 2.0*randgen() - 1.0;
		double w = randnum1*randnum1 + randnum2*randnum2;
		while(w >= 1.0){ 
			randnum1 = 2.0*randgen() - 1.0;
			randnum2 = 2.0*randgen() - 1.0;
			w = randnum1*randnum1 + randnum2*randnum2;
		}
		w = pow(-double(2.0)*log(w)/w,double(0.5));
		cached_value = randnum2*w;
		value_is_cached = 1;
		return(randnum1*w);
	}
	cout << "gaussrandgen, no cached choice error" << endl;
	return(0.0);
}
