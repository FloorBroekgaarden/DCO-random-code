#include "lambdaNanjing.h"

double lambdaNanjing(double Mzamsa, double Ma, double Mca, double Ra, int stellarType, double metallicity, const programOptions &options){
	/*
	Function for calculating the common envelope lambda parameter using the "Nanjing" prescription
	from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) as implemented in STARTRACK

	This function is copied from STARTRACK courtesy of Chris Belczynski.

	Parameters
	-----------
	Mzamsa : double
		Zero age main sequence mass of star in solar masses
	Ma : double
		Current mass of star in solar masses
	Ra : double
		Current radius of star in solar radii
	stellarType : int
		Current stellar type (according to Hurley et al 2000)
	metallicity : double
		Metallicity of star
	options : programOptions
		User specified program options

	Options
	--------
	lambda : double
		Nanjing lambda for use in common envelope
	*/
	bool	debugging = false;
	double lambda,lamb,lamg,maxb,maxg;
	double a0,a1,a2,a3,a4,a5;
	double b0,b1,b2,b3,b4,b5;
	double x,y1,y2,Zlimit,Menva;
	double Rmin,Rmax;
	int manual;

	lambda		= 1.0;					// Default in STARTRACK
	double ZZsun		= 0.02;					// Default in STARTRACK
	double	ZZ 			= metallicity;

    if(debugging){
        std::cout << "Inside lambdaNanjing function" << std::endl;
        std::cout << "\nLambda:\t" << lambda << std::endl;
    }

	if(stellarType == NAKED_HELIUM_STAR_HETZSPRUNG_GAP or stellarType == NAKED_HELIUM_STAR_GIANT_BRANCH) {
		/* For Helium stars: always use Natasha's fit */

		if(debugging){std::cout << "\nHelium stars -- Natasha fit\n" << std::endl;}

		Rmin = 0.25;			/* minimum considered radius: Natasha */
		Rmax = 120.0;			/* maximum considered radius: Natasha */
	
		if(Ra < Rmin){
	 		lambda=0.3*pow(Rmin,-0.8);
	 	}
		else if(Ra > Rmax){
	 		lambda=0.3*pow(Rmax,-0.8);
	 	}
		else{
	 		lambda=0.3*pow(Ra,-0.8);
	 	}
	}
	else if(stellarType >= HERTZSPRUNG_GAP and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH){

        // Use Nanjing for H-rich stars
        Zlimit = 0.0105;
        Menva  = Ma - Mca;					/* current He core mass . SIMON: Seems more like envelope mass to me ...*/
        manual = 0;                        	/* then approximated by hand: manual=1 */

        if(ZZ > Zlimit){                                                            /* Z>0.5 Zsun: popI */
            if(Mzamsa < 1.5){

                maxb = 2.5;
                maxg = 1.5;

                if(Ra > 200.0) {
                    lamb   = 0.05;
                    lamg   = 0.05;
                    manual = 1;
                }
                else if((stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) and Ra > 2.7){
                    lamb=-Ra*9.18e-03+2.33;
                    lamg=-Ra*4.59e-03+1.12;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH){
                    a0=8.35897;  a1=-18.89048; a2=10.47651; a3=0.99352;  a4=0.0; a5=0.0;  /* coeffs: with internal */
                    b0=17.58328; b1=-34.84355; b2=10.70536; b3=8.49042;  b4=0.0; b5=0.0;  /* coeffs: no internal */
                }
                else if(stellarType == CORE_HELIUM_BURNING){
                    a0=46.00978; a1=-298.64993; a2=727.40936; a3=-607.66797; a4=0.0; a5=0.0;
                    b0=63.61259; b1=-399.89494; b2=959.62055; b3=-795.20699; b4=0.0; b5=0.0;
                }
                else {
                    lamb=-Ra*3.57e-04+0.1;
                    lamg=-Ra*3.57e-04+0.1;
                    manual=1;
                }
            }
            else if(Mzamsa<2.5) {
           
                maxb=4.0;
                maxg=2.0;

                if(Ra>340.0) {
                    lamb=3.589970;
                    lamg=0.514132;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra > 8.5 and Ra < 60.0) {
                    lamb=3.0;
                    lamg=1.2;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=2.05363; a1=-0.00685; a2=-3.42739e-04; a3=3.93987e-06; a4=-1.18237e-08; a5=0.0;
                    b0=1.07658; b1=-0.01041; b2=-4.90553e-05; b3=1.13528e-06; b4=-3.91609e-09; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=34.41826; a1=-6.65259; a2=0.43823; a3=-0.00953; a4=0.0; a5=0.0;
                    b0=13.66058; b1=-2.48031; b2=0.15275; b3=-0.00303; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.88954; a1=0.0098;  a2=-3.1411e-05;  a3=7.66979e-08; a4=0.0;         a5=0.0;
                    b0=0.48271; b1=0.00584; b2=-6.22051e-05; b3=2.41531e-07; b4=-3.1872e-10; b5=0.0;
                }
            }
            else if(Mzamsa<3.5) {

                maxb=500.0;
                maxg=10.0;

                if(Ra>400.0) {
                    lamb=116.935557;
                    lamg=0.848808;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=2.40831; a1=-0.42459; a2=0.03431; a3=-9.26879e-04; a4=8.24522e-06; a5=0.0;
                    b0=1.30705; b1=-0.22924; b2=0.01847; b3=-5.06216e-04; b4=4.57098e-06; b5=0.0;
                    maxb=2.5; maxg=1.5;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-42.98513; a1=7.90134; a2=-0.54646; a3=0.01863; a4=-3.13101e-04; a5=2.07468e-06;
                    b0=-6.73842;  b1=1.06656; b2=-0.05344; b3=0.00116; b4=-9.34446e-06; b5=0.0;
                    maxb=2.5; maxg=1.5;
                }
                else {
                    a0=-0.04669; a1=0.00764; a2=-4.32726e-05; a3=9.31942e-08; a4=0.0;         a5=0.0;
                    b0=0.44889;  b1=0.01102; b2=-6.46629e-05; b3=5.66857e-09; b4=7.21818e-10; b5=-1.2201e-12;
                }
            }
            else if(Mzamsa<4.5) {

                maxb=1000.0;
                maxg=8.0;

                if(Ra>410.0) {
                    lamb=52.980056;
                    lamg=1.109736;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.8186;  a1=-0.17464; a2=0.00828; a3=-1.31727e-04; a4=7.08329e-07; a5=0.0;
                    b0=1.02183; b1=-0.1024;  b2=0.00493; b3=-8.16343e-05; b4=4.55426e-07; b5=0.0;
                    maxb=2.5; maxg=1.5;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-7.3098;  a1=0.56647; a2=-0.01176; a3=7.90112e-05; a4=0.0; a5=0.0;
                    b0=-3.80455; b1=0.29308; b2=-0.00603; b3=4.00471e-05; b4=0.0; b5=0.0;
                    maxb=2.5; maxg=1.5;
                }
                else {
                    a0=-0.37322; a1=0.00943; a2=-3.26033e-05; a3=5.37823e-08; a4=0.0; a5=0.0;
                    b0=0.13153;  b1=0.00984; b2=-2.89832e-05; b3=2.63519e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<5.5) {

                maxb=1000.0;
                maxg=8.0;

                if(Ra>430.0) {
                    lamb=109.593522;
                    lamg=1.324248;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.52581; a1=-0.08125; a2=0.00219; a3=-2.0527e-05;  a4=6.79169e-08; a5=0.0;
                    b0=0.85723; b1=-0.04922; b2=0.00137; b3=-1.36163e-05; b4=4.68683e-08; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-9.93647; a1=0.42831; a2=-0.00544; a3=2.25848e-05; a4=0.0; a5=0.0;
                    b0=-5.33279; b1=0.22728; b2=-0.00285; b3=1.16408e-05; b4=0.0; b5=0.0;
                }
                else {
                    a0=-0.80011; a1=0.00992; a2=-3.03247e-05; a3=5.26235e-08;  a4=0.0; a5=0.0;
                    b0=-0.00456; b1=0.00426; b2=4.71117e-06;  b3=-1.72858e-08; b4=0.0; b5=0.0;
                }                
            }
            else if(Mzamsa<6.5) {

                maxb=25.5;
                maxg=5.0;

                if(Ra>440.0) {
                    lamb=16.279603;
                    lamg=1.352166;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.41601; a1=-0.04965; a2=8.51527e-04; a3=-5.54384e-06; a4=1.32336e-08; a5=0.0;
                    b0=0.78428; b1=-0.02959; b2=5.2013e-04;  b3=-3.45172e-06; b4=8.17248e-09; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=13.91465; a1=-0.55579; a2=0.00809; a3=-4.94872e-05; a4=1.08899e-07; a5=0.0;
                    b0=7.68768;  b1=-0.30723; b2=0.00445; b3=-2.70449e-05; b4=5.89712e-08; b5=0.0;
                }
                else {
                    a0=-2.7714; a1=0.06467;  a2=-4.01537e-04; a3=7.98466e-07;  a4=0.0; a5=0.0;
                    b0=0.23083; b1=-0.00266; b2=2.21788e-05;  b3=-2.35696e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<7.5) {

                maxb=9.0;
                maxg=3.0;

                if(Ra>420.0) {
                    lamb=5.133959;
                    lamg=1.004036;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.38344; a1=-0.04093; a2=5.78952e-04; a3=-3.19227e-06; a4=6.40902e-09; a5=0.0;
                    b0=0.76009; b1=-0.02412; b2=3.47104e-04; b3=-1.92347e-06; b4=3.79609e-09; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=4.12387; a1=-0.12979; a2=0.00153;     a3=-7.43227e-06; a4=1.29418e-08; a5=0.0;
                    b0=2.18952; b1=-0.06892; b2=8.00936e-04; b3=-3.78092e-06; b4=6.3482e-09;  b5=0.0;
                }
                else {
                    a0=-0.63266; a1=0.02054;  a2=-1.3646e-04; a3=2.8661e-07;   a4=0.0; a5=0.0;
                    b0=0.26294;  b1=-0.00253; b2=1.32272e-05; b3=-7.12205e-09; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<8.5) {

                maxb=7.0;
                maxg=3.0;

                if(Ra>490.0) {
                    lamb=4.342985;
                    lamg=0.934659;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.35516; a1=-0.03414; a2=4.02065e-04; a3=-1.85931e-06; a4=3.08832e-09; a5=0.0;
                    b0=0.73826; b1=-0.01995; b2=2.37842e-04; b3=-1.09803e-06; b4=1.79044e-09; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-3.89189; a1=0.19378; a2=-0.0032;  a3=2.39504e-05; a4=-8.28959e-08; a5=1.07843e-10;
                    b0=-2.24354; b1=0.10918; b2=-0.00179; b3=1.33244e-05; b4=-4.57829e-08; b5=5.90313e-11;
                    maxb=1.0; maxg=0.5;
                }
                else {
                    a0=-0.1288; a1=0.0099;   a2=-6.71455e-05;  a3=1.33568e-07;   a4=0.0; a5=0.0;
                    b0=0.26956; b1=-0.00219; b2=7.97743e-06;   b3=-1.53296e-09;  b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<9.5) {

                maxb=4.0;
                maxg=2.0;

                if(Ra>530.0) {
                    lamb=2.441672;
                    lamg=0.702310;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.32549; a1=-0.02845; a2=2.79097e-04; a3=-1.07254e-06; a4=1.46801e-09; a5=0.0;
                    b0=0.71571; b1=-0.01657; b2=1.64607e-04; b3=-6.31935e-07; b4=8.52082e-10; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.86369; a1=-0.00995; a2=4.80837e-05;  a3=-6.10454e-08; a4=-2.79504e-12; a5=0.0;
                    b0=-0.7299; b1=0.0391;   b2=-5.78132e-04; b3=3.7072e-06;   b4=-1.07036e-08; b5=1.14833e-11;
                }
                else {
                    a0=1.19804; a1=-0.01961; a2=1.28222e-04; a3=-3.41278e-07; a4=3.35614e-10; a5=0.0;
                    b0=0.40587; b1=-0.0051;  b2=2.73866e-05; b3=-5.74476e-08; b4=4.90218e-11; b5=0.0;
                }
            }
            else if(Mzamsa<11.0) {

                maxb=3.0;
                maxg=1.5;

                if(Ra>600.0) {
                    lamb=1.842314;
                    lamg=0.593854;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.29312; a1=-0.02371; a2=1.93764e-04; a3=-6.19576e-07; a4=7.04227e-10; a5=0.0;
                    b0=0.69245; b1=-0.01398; b2=1.17256e-04; b3=-3.81487e-07; b4=4.35818e-10; b5=0.0;
                    maxb=1.0; maxg=0.6;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.74233; a1=-0.00623; a2=2.04197e-05; a3=-1.30388e-08; a4=0.0; a5=0.0;
                    b0=0.36742; b1=-0.00344; b2=1.27838e-05; b3=-1.0722e-08;  b4=0.0; b5=0.0;
                }
                else {
                    a0=0.3707;  a1=2.67221e-04; a2=-9.86464e-06; a3=2.26185e-08; a4=0.0; a5=0.0;
                    b0=0.25549; b1=-0.00152;    b2=3.35239e-06;  b3=2.24224e-10; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<13.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>850.0) {
                    lamb=0.392470;
                    lamg=0.176660;
                    manual=1;
                }
                else if(Ra>0.0 and Ra<=350.0) {
                    a0=1.28593; a1=-0.02209; a2=1.79764e-04; a3=-6.21556e-07; a4=7.59444e-10; a5=0.0;
                    b0=0.68544; b1=-0.01394; b2=1.20845e-04; b3=-4.29071e-07; b4=5.29169e-10; b5=0.0;
                }
                else if(Ra>350.0 and Ra<=600.0) {
                    a0=-11.99537; a1=0.0992;  a2=-2.8981e-04; a3=3.62751e-07;  a4=-1.65585e-10; a5=0.0;
                    b0=0.46156;   b1=-0.0066; b2=3.9625e-05;  b3=-9.98667e-08; b4=-8.84134e-11; b5=0.0;
                }
                else {
                    a0=-58.03732; a1=0.23633; a2=-3.20535e-04; a3=1.45129e-07; a4=0.0; a5=0.0;
                    b0=-15.11672; b1=0.06331; b2=-8.81542e-05; b3=4.0982e-08;  b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<15.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>1000.0) {
                    lamb=0.414200;
                    lamg=0.189008;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra > 69.0 and Ra < 126.0) {
                    lamb=Ra*(-8.77e-04)+0.5;
                    lamg=0.18;
                    manual=1;
                }
                else if((stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) and Ra>190.0 and Ra<600.0) {
                    lamb=0.15;
                    lamg=0.15;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.39332;  a1=-0.0318; a2=3.95917e-04; a3=-2.23132e-06; a4=4.50831e-09; a5=0.0;
                    b0=0.78215; b1=-0.02326; b2=3.25984e-04; b3=-1.94991e-06; b4=4.08044e-09; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=1.12889; a1=-0.00901; a2=3.04077e-05; a3=-4.31964e-08; a4=2.14545e-11; a5=0.0;
                    b0=0.568;   b1=-0.0047;  b2=1.57818e-05; b3=-2.21207e-08; b4=1.08472e-11; b5=0.0;
                }
                else {
                    a0=-106.90553; a1=0.36469; a2=-4.1472e-04;  a3=1.57349e-07; a4=0.0; a5=0.0;
                    b0=-39.93089;  b1=0.13667; b2=-1.55958e-04; b3=5.94076e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<18.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>1050.0) {
                    lamb=0.2;
                    lamg=0.1;
                    manual=1;
                }
                else if((stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) and Ra>120.0 and Ra<170.0) {
                    lamb=0.2;
                    lamg=0.2;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.43177; a1=-0.03533; a2=5.11128e-04; a3=-3.57633e-06; a4=9.36778e-09; a5=0.0;
                    b0=0.85384; b1=-0.03086; b2=5.50878e-04; b3=-4.37671e-06; b4=1.25075e-08; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.84143; a1=-0.00576; a2=1.68854e-05; a3=-2.0827e-08;  a4=8.97813e-12; a5=0.0;
                    b0=0.36014; b1=-0.00254; b2=7.49639e-06; b3=-9.20103e-09; b4=3.93828e-12; b5=0.0;
                }
                else {
                    a0=-154.70559; a1=0.46718; a2=-4.70169e-04; a3=1.57773e-07; a4=0.0; a5=0.0;
                    b0=-65.39602;  b1=0.19763; b2=-1.99078e-04; b3=6.68766e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<35.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>1200.0) {
                    lamb=0.05;
                    lamg=0.05;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    lamb=1.2*exp(-Ra/90.0);
                    lamg=0.55*exp(-Ra/160.0);
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.48724; a1=-0.00177;    a2=2.60254e-06; a3=-1.25824e-09; a4=0.0; a5=0.0;
                    b0=0.22693; b1=-8.7678e-04; b2=1.28852e-06; b3=-6.12912e-10; b4=0.0; b5=0.0;
                }
                else {
                    a0=-260484.85724; a1=4.26759e+06; a2=-2.33016e+07; a3=4.24102e+07; a4=0.0; a5=0.0;
                    b0=-480055.67991; b1=7.87484e+6;  b2=-4.30546e+07; b3=7.84699e+07; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<75.0) {
                maxb=1.0;
                maxg=0.5;

               a0=0.31321; a1=-7.50384e-04; a2=5.38545e-07; a3=-1.16946e-10; a4=0.0; a5=0.0;
               b0=0.159;   b1=-3.94451e-04; b2=2.88452e-07; b3=-6.35132e-11; b4=0.0; b5=0.0;
            }
            else {
                maxb=1.0;
                maxg=0.5;

                a0=0.376;  a1=-0.0018;  a2=2.81083e-06; a3=-1.67386e-09; a4=3.35056e-13; a5=0.0;
                b0=0.2466; b1=-0.00121; b2=1.89029e-06; b3=-1.12066e-09; b4=2.2258e-13;  b5=0.0;
            }
        }
        else {                                  /* Z<=0.5 Zsun: popI and popII */
            if(Mzamsa<1.5) {
                maxb=2.0;
                maxg=1.5;

                if(Ra>160.0) {
                    lamb=0.05;
                    lamg=0.05;
                    manual=1;
                }
                else if((stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) and Ra>12.0) {
                    lamb=1.8*exp(-Ra/80.0);
                    lamg=exp(-Ra/45.0);
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.37294; a1=-0.05825; a2=0.00375; a3=-7.59191e-05; a4=0.0; a5=0.0;
                    b0=0.24816; b1=-0.04102; b2=0.0028;  b3=-6.20419e-05; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.24012; a1=-0.01907; a2=6.09529e-04; a3=-8.17819e-06; a4=4.83789e-08; a5=-1.04568e-10;
                    b0=0.15504; b1=-0.01238; b2=3.96633e-04; b3=-5.3329e-06;  b4=3.16052e-08; b5=-6.84288e-11;
                }
            }
            else if(Mzamsa<2.5) {
                maxb=4.0;
                maxg=2.0;

                if(Ra>350.0) {
                    lamb=2.868539;
                    lamg=0.389991;
                    manual=1;
                }
                else if((stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) and Ra>22.0 and Ra<87.0) {
                    lamb=1.95;
                    lamg=0.85;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>6.0 and Ra<50.0) {
                    lamb=0.8;
                    lamg=0.35;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=2.56108; a1=-0.75562; a2=0.1027;  a3=-0.00495; a4=8.05436e-05; a5=0.0;
                    b0=1.41896; b1=-0.4266;  b2=0.05792; b3=-0.00281; b4=4.61e-05;    b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-103.92538; a1=25.37325; a2=-2.03273; a3=0.0543;  a4=0.0; a5=0.0;
                    b0=-56.03478;  b1=13.6749;  b2=-1.09533; b3=0.02925; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.5452;  a1=0.00212;      a2=6.42941e-05; a3=-1.46783e-07; a4=0.0;        a5=0.0;
                    b0=0.30594; b1=-9.58858e-04; b2=1.12174e-04; b3=-1.04079e-06; b4=3.4564e-09; b5=-3.91536e-12;
                }
            }
            else if(Mzamsa<3.5) {
                maxb=600.0;
                maxg=2.0;

                if(Ra>400.0) {
                    lamb=398.126442;
                    lamg=0.648560;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>36.0 and Ra<53.0) {
                    lamb=1.0;
                    lamg=1.0;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.7814;  a1=-0.17138; a2=0.00754; a3=-9.02652e-05; a4=0.0; a5=0.0;
                    b0=0.99218; b1=-0.10082; b2=0.00451; b3=-5.53632e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-12.40832; a1=1.59021; a2=-0.06494; a3=8.69587e-04; a4=0.0; a5=0.0;
                    b0=-6.47476;  b1=0.8328;  b2=-0.03412; b3=4.58399e-04; b4=0.0; b5=0.0;
                }
                else {
                    a0=-0.475;  a1=-0.00328; a2=1.31101e-04; a3=-6.03669e-07; a4=8.49549e-10; a5=0.0;
                    b0=0.05434; b1=0.0039;   b2=9.44609e-06; b3=-3.87278e-08; b4=0.0;         b5=0.0;
                }
            }
            else if(Mzamsa<4.5) {
                maxb=600.0;
                maxg=2.0;

                if(Ra>410.0) {
                    lamb=91.579093;
                    lamg=1.032432;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>19.0 and Ra<85.0) {
                    lamb=0.255;
                    lamg=0.115;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.65914; a1=-0.10398; a2=0.0029;  a3=-2.24862e-05; a4=0.0; a5=0.0;
                    b0=0.92172; b1=-0.06187; b2=0.00177; b3=-1.42677e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-5.89253; a1=0.54296; a2=-0.01527; a3=1.38354e-04; a4=0.0; a5=0.0;
                    b0=-3.21299; b1=0.29583; b2=-0.00833; b3=7.55646e-05; b4=0.0; b5=0.0;
                }
                else {
                    a0=-0.2106; a1=-0.01574; a2=2.01107e-04; a3=-6.90334e-07; a4=7.92713e-10; a5=0.0;
                    b0=0.36779; b1=-0.00991; b2=1.19411e-04; b3=-3.59574e-07; b4=3.33957e-10; b5=0.0;
                }
            }
            else if(Mzamsa<5.5) {
                maxb=10.0;
                maxg=3.0;

                if(Ra>320) {
                    lamb=7.618019;
                    lamg=1.257919;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>85.0 and Ra<120.0) {
                    lamb=0.4;
                    lamg=0.1;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.58701; a1=-0.06897; a2=0.00129;     a3=-6.99399e-06; a4=0.0; a5=0.0;
                    b0=0.87647; b1=-0.04103; b2=7.91444e-04; b3=-4.41644e-06; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=-0.67176; a1=0.07708; a2=-0.00175;    a3=1.1991e-05;  a4=0.0; a5=0.0;
                    b0=-0.38561; b1=0.0427;  b2=-9.6948e-04; b3=6.64455e-06; b4=0.0; b5=0.0;
                }
                else {
                    a0=-0.12027; a1=0.01981;  a2=-2.27908e-04; a3=7.55556e-07;  a4=0.0; a5=0.0;
                    b0=0.31252;  b1=-0.00527; b2=3.60348e-05;  b3=-3.22445e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<6.5) {
                maxb=4.0;
                maxg=1.5;

                if(Ra>330.0) {
                    lamb=2.390575;
                    lamg=0.772091;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>115.0 and Ra<165.0) {
                    lamb=0.2;
                    lamg=0.1;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.527;   a1=-0.04738; a2=6.1373e-04;  a3=-2.36835e-06; a4=0.0; a5=0.0;
                    b0=0.83636; b1=-0.02806; b2=3.73346e-04; b3=-1.47016e-06; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.30941; a1=0.00965; a2=-2.31975e-04; a3=1.26273e-06; a4=0.0; a5=0.0;
                    b0=0.14576; b1=0.00562; b2=-1.30273e-04; b3=7.06459e-07; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.26578; a1=0.00494;  a2=-7.02203e-05; a3=2.25289e-07; a4=0.0; a5=0.0;
                    b0=0.26802; b1=-0.00248; b2=6.45229e-06;  b3=1.69609e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<7.5) {
                maxb=2.5;
                maxg=1.0;

                if(Ra>360.0) {
                    lamb=1.878174;
                    lamg=0.646353;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>150.0 and Ra<210.0) {
                    lamb=0.2;
                    lamg=0.1;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.49995; a1=-0.03921; a2=4.2327e-04; a3=-1.37646e-06; a4=0.0; a5=0.0;
                    b0=0.81688; b1=-0.02324; b2=2.5804e-04; b3=-8.54696e-07; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.44862; a1=0.00234; a2=-9.23152e-05; a3=4.67797e-07; a4=0.0; a5=0.0;
                    b0=0.21873; b1=0.00154; b2=-5.18806e-05; b3=2.60283e-07; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.8158;  a1=-0.01633; a2=1.46552e-04; a3=-5.75308e-07; a4=8.77711e-10; a5=0.0;
                    b0=0.26883; b1=-0.00219; b2=4.12941e-06; b3=1.33138e-08;  b4=0.0;         b5=0.0;
                }
            }
            else if(Mzamsa<8.5) {
                maxb=2.0;
                maxg=1.0;

                if(Ra>400.0) {
                    lamb=1.517662;
                    lamg=0.553169;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>190.0 and Ra<260.0) {
                    lamb=0.2;
                    lamg=0.1;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.46826; a1=-0.03184; a2=2.85622e-04; a3=-7.91228e-07; a4=0.0; a5=0.0;
                    b0=0.79396; b1=-0.01903; b2=1.77574e-04; b3=-5.04262e-07; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.50221; a1=-3.19021e-04; a2=-3.81717e-05; a3=1.80726e-07; a4=0.0; a5=0.0;
                    b0=0.24748; b1=-9.9338e-05;  b2=-1.99272e-05; b3=9.47504e-08; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.74924; a1=-0.01233; a2=9.55715e-05; a3=-3.37117e-07; a4=4.67367e-10; a5=0.0;
                    b0=0.25249; b1=-0.00161; b2=8.35478e-07; b3=1.25999e-08;  b4=0.0;         b5=0.0;
                }
            }
            else if(Mzamsa<9.5) {
                maxb=1.6;
                maxg=1.0;

                if(Ra>440.0) {
                    lamb=1.136394;
                    lamg=0.478963;
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING and Ra>180.0 and Ra<300.0) {
                    lamb=0.15;
                    lamg=0.08;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.49196; a1=-0.03247; a2=3.08066e-04; a3=-9.53247e-07; a4=0.0; a5=0.0;
                    b0=0.805;   b1=-0.02;    b2=2.01872e-04; b3=-6.4295e-07;  b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.39342; a1=0.00259;     a2=-4.97778e-05; a3=1.69533e-07; a4=0.0; a5=0.0;
                    b0=0.20796; b1=6.62921e-04; b2=-1.84663e-05; b3=6.58983e-08; b4=0.0; b5=0.0;
                }
                else {
                    a0=0.73147; a1=-0.01076; a2=7.54308e-05; a3=-2.4114e-07;  a4=2.95543e-10; a5=0.0;
                    b0=0.31951; b1=-0.00392; b2=2.31815e-05; b3=-6.59418e-08; b4=7.99575e-11; b5=0.0;
                }
            }
            else if(Mzamsa<11.0) {
                maxb=1.6;
                maxg=1.0;

                if(Ra>500.0) {
                    lamb=1.068300;
                    lamg=0.424706;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    lamb=1.75*exp(-Ra/35.);
                    lamg=0.9*exp(-Ra/35.);
                    manual=1;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.75746; a1=-0.00852; a2=3.51646e-05; a3=-4.57725e-08; a4=0.0; a5=0.0;
                    b0=0.35355; b1=-0.00388; b2=1.56573e-05; b3=-1.98173e-08; b4=0.0; b5=0.0;
                }
                else {
                    a0=-9.26519; a1=0.08064;  a2=-2.30952e-04; a3=2.21986e-07; a4=0.0; a5=0.0;
                    b0=0.81491;  b1=-0.00161; b2=-8.13352e-06; b3=1.95775e-08; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<13.0) {
                maxb=1.6;
                maxg=1.0;

                if(Ra>600.0) {
                    lamb=0.537155;
                    lamg=0.211105;
                    manual=1;
                }
                else if((stellarType == CORE_HELIUM_BURNING and Ra>200.0 and Ra<410.0) or (stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>390.0 and Ra<460.0)) {
                    lamb=0.08;
                    lamg=0.05;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.63634; a1=-0.04646; a2=7.49351e-04; a3=-5.23622e-06; a4=0.0; a5=0.0;
                    b0=1.17934; b1=-0.08481; b2=0.00329;     b3=-4.69096e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.85249; a1=-0.00861; a2=2.99246e-05; a3=-3.21416e-08; a4=0.0; a5=0.0;
                    b0=0.37188; b1=-0.00365; b2=1.24944e-05; b3=-1.32388e-08; b4=0.0; b5=0.0;
                }
                else {
                    a0=-51.15252; a1=0.30238; a2=-5.95397e-04; a3=3.91798e-07; a4=0.0; a5=0.0;
                    b0=-13.44;    b1=0.08141; b2=-1.641e-04;   b3=1.106e-07;   b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<15.0) {
                maxb=1.6;
                maxg=1.0;

                if(Ra>650.0) {
                    lamb=0.3;
                    lamg=0.160696;
                    manual=1;
                }
                else if((stellarType == CORE_HELIUM_BURNING and Ra>250.0 and Ra<490.0) or (stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>480.0 and Ra<540.0)) {
                    lamb=0.06;
                    lamg=0.05;
                    manual=1;
                }
                else if(stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>540.0 and Ra<650.0) {
                    lamb=Ra*1.8e-03-0.88;
                    lamg=Ra*9.1e-04-0.43;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.45573; a1=-0.00937; a2=-0.00131; a3=3.07004e-05;  a4=0.0; a5=0.0;
                    b0=1.19526; b1=-0.08503; b2=0.00324;  b3=-4.58919e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.85271; a1=-0.00793; a2=2.5174e-05;  a3=-2.4456e-08;  a4=0.0; a5=0.0;
                    b0=0.36163; b1=-0.00328; b2=1.03119e-05; b3=-9.92712e-09; b4=0.0; b5=0.0;
                }
                else {
                    a0=-140.0;   a1=0.7126;  a2=-0.00121;     a3=6.846e-07;   a4=0.0; a5=0.0;
                    b0=-44.1964; b1=0.22592; b2=-3.85124e-04; b3=2.19324e-07; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<18.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>750.0) {
                    lamb=0.5;
                    lamg=0.204092;
                    manual=1;
                }
                else if((stellarType == CORE_HELIUM_BURNING and Ra>200.0 and Ra<570.0) or (stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>560.0 and Ra<650.0)) {
                    lamb=0.1;
                    lamg=0.05;
                    manual=1;
                }
                else if(stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>650.0 and Ra<750.0) {
                    lamb=Ra*4.0e-03-2.5;
                    lamg=Ra*1.5e-03-0.93;
                    manual=1;
                }
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.33378; a1=0.01274;  a2=-0.00234; a3=4.6036e-05;   a4=0.0; a5=0.0;
                    b0=1.17731; b1=-0.07834; b2=0.00275;  b3=-3.58108e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.83254; a1=-0.00696; a2=1.9597e-05;  a3=-1.67985e-08; a4=0.0; a5=0.0;
                    b0=0.34196; b1=-0.0028;  b2=7.82865e-06; b3=-6.66684e-09; b4=0.0; b5=0.0;
                }
                else {
                    a0=-358.4;     a1=1.599;   a2=-0.00238;    a3=1.178e-06;   a4=0.0; a5=0.0;
                    b0=-118.13757; b1=0.52737; b2=-7.8479e-04; b3=3.89585e-07; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<35.0) {
                maxb=1.5;
                maxg=1.0;

                if(Ra>900.0) {
                    lamb=0.2;
                    lamg=0.107914;
                    manual=1;
                }
                else if((stellarType == CORE_HELIUM_BURNING and Ra>230.0 and Ra<755.0) or (stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>725.0 and Ra<850.0)) {
                    lamb=0.1;
                    lamg=0.05;
                    manual=1;
                }
                else if(stellarType == EARLY_ASYMPTOTIC_GIANT_BRANCH and Ra>850.0 and Ra<900.0) {
                    lamb=Ra*2.0e-03-1.6;
                    lamg=Ra*1.0e-03-0.8;
                    manual=1;
                } 
                else if(stellarType == HERTZSPRUNG_GAP or stellarType == FIRST_GIANT_BRANCH) {
                    a0=1.27138; a1=0.00538;  a2=-0.0012; a3=1.80776e-05;  a4=0.0; a5=0.0;
                    b0=1.07496; b1=-0.05737; b2=0.00153; b3=-1.49005e-05; b4=0.0; b5=0.0;
                }
                else if(stellarType == CORE_HELIUM_BURNING) {
                    a0=0.69746; a1=-0.0043;  a2=8.97312e-06; a3=-5.83402e-09; a4=0.0; a5=0.0;
                    b0=0.26691; b1=-0.00161; b2=3.3378e-06;  b3=-2.1555e-09;  b4=0.0; b5=0.0;
                }
                else {
                    a0=-436.00777; a1=1.41375; a2=-0.00153;     a3=5.47573e-07; a4=0.0; a5=0.0;
                    b0=-144.53456; b1=0.46579; b2=-4.99197e-04; b3=1.78027e-07; b4=0.0; b5=0.0;
                }
            }
            else if(Mzamsa<75.0) {
                maxb=20.0;
                maxg=3.0;

                a0=0.821;   a1=-0.00669; a2=1.57665e-05; a3=-1.3427e-08;  a4=3.74204e-12; a5=0.0;
                b0=0.49287; b1=-0.00439; b2=1.06766e-05; b3=-9.22015e-09; b4=2.58926e-12; b5=0.0;
            }
            else {
                maxb=4.0;
                maxg=2.0;

                a0=1.25332; a1=-0.02065; a2=1.3107e-04;  a3=-3.67006e-07; a4=4.58792e-10; a5=-2.09069e-13;
                b0=0.81716; b1=-0.01436; b2=9.31143e-05; b3=-2.6539e-07;  b4=3.30773e-10; b5=-1.51207e-13;
            }
        }

        if(manual==0) {
            if(ZZ>Zlimit and (Mzamsa<1.5 or (Mzamsa>=18.0 and Mzamsa<25.0 and stellarType >= EARLY_ASYMPTOTIC_GIANT_BRANCH and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH))) {
                x=Menva/Ma;
                y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
                y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
                lamb=1.0/y1;
                lamg=1.0/y2;
            }
            else if((ZZ>Zlimit and Mzamsa>=2.5 and Mzamsa<5.5 and stellarType >= EARLY_ASYMPTOTIC_GIANT_BRANCH and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH) or (ZZ<=Zlimit and Mzamsa>=2.5 and Mzamsa<4.5 and stellarType >= EARLY_ASYMPTOTIC_GIANT_BRANCH and stellarType <= THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH)) {
                x=Ra;
                y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
                y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
                lamb=pow(10.0,y1);
                lamg=y2;
            }
            else {
                x=Ra;
                y1=a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x+a5*x*x*x*x*x;
                y2=b0+b1*x+b2*x*x+b3*x*x*x+b4*x*x*x*x+b5*x*x*x*x*x;
                lamb=y1;
                lamg=y2;
            }
        }

        // Limit lambda to some 'reasonable' range
        if(lamb<0.05){
            lamb=0.05;
        }
        if(lamg<0.05){
            lamg=0.05;
        }
        if(lamb>maxb){
            lamb=maxb;
        }
        if(lamg>maxg){
            lamg=maxg;
        }

        // Calculate lambda as some combination of lambda_b and lambda_g by
        // lambda = alpha_th • lambda_b    +  (1-alpha_th) • lambda_g 
        // Note that this is different to STARTRACK
        lambda = (options.commonEnvelopeAlphaThermal * lamb) + (1.0 - options.commonEnvelopeAlphaThermal)*lamg;

            //GOT TO HERE
    }
    else{
        // Cannot compute lambdaNanjing for stellar types 0,1,7,11,12,13,14
        if(options.debugging){
            std::cerr << "Error: calling lambdaNanjing for a star with stellar type " << stellarType << " not supported." << std::endl;
        }
    }
    
    if(debugging)
        std::cout << "\nLambda:\t" << lambda << std::endl;

    return lambda;

}