#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NR_RNDNUM 365		//random generator: delicate mode
                                // undetectable correlations
int icall;                      // used by random number generator	
int	iv[NR_RNDNUM];		// used by random number generator


#define quick_and_dirty_rand() \
	rand()/double(RAND_MAX)


void rini(int iran);		// random generator initialisation
double rndnum();		// random generator: uniformly, between 0 and 1.
//Simple template of a main function using rndnum()
//int main(){
//	int i;
//      unsigned int rseed;
//
//	// initialisation ("set seed"): 0<= rseed < 259200
//	rini( rseed );
//
//      // a random number is a double returned by rndnum()
//	for (i=0;i<1000;i++)
//		printf("%.7f\t",rndnum());
//	printf("\n");
//
//	return 0;
//}



//
//  file random.f, Ulli Wolff, 15.9.95
//
//  random generator modified RCARRY, 
//  see M. Luescher, Comp. Phys. Comm. 79(1994)100
//  and literature quoted there
//
//  This generator is relatively slow, but of KNOWN good quality
//  The speed strongly depends on the value ithrow
//  One way to speed up would be vectoriztion over parallel copies
//
//  ALWAYS CHANGE PARAMETERS IN BOTH ROUTINES OF THIS FILE
//
//  values for ithrow: 0 for simple uses, 24 for normal, 199 for
//  practically no detectable correlations, 365 for tests in delicated cases;
//  the bigger ithrow, the slower the generator!
//
//  0 <= rndnum < 1  (rndnum=0 DOES occur occasionally)
//  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
double rndnum(){
		


      int ir,is,irs,ir1,ikeep,ithrow,ivl1,ivl;
	  int i,j,ivn;

	  long ibase;
      double basein;


	  double frndnum;

		ir=24;
		is=10;
		irs=ir-is;
		ir1=ir-1;
		ikeep=ir;

		ithrow=24;			// normal

		ivl1=ikeep+ithrow;
		ivl=ivl1+ir1;		// this is the size of iv: ithrow + 2*ir - 1 (+1)
		ibase=long(pow(2,24));
		basein=1.0/ibase;


//  state of the generator (to be saved and reloaded for continuation):
//		common/rando/ iv(0:ivl),icall
//      save /rando/

	if (icall==ikeep){
//  disregard ithrow numbers:
        for(j=1;j<=ithrow;j++){
			ivn=iv[irs+icall]-iv[icall];
          icall=icall+1;
//  carry bit:
			if (ivn<0){
				iv[icall]=iv[icall]+1;
				ivn=ivn+ibase;
			}

			iv[icall+ir1]=ivn;
		}
//  copy last ir numbers to the beginning:
        for(i=0;i<=ir1;i++)
	          iv[i]=iv[ivl1+i];
		
		icall=0;
    }

//  basic cycle:
      ivn=iv[irs+icall]-iv[icall];
      icall=icall+1;
//  carry bit:
      if (ivn<0){
        iv[icall]=iv[icall]+1;
        ivn=ivn+ibase;
      }

      iv[icall+ir1]=ivn;

//  convert to floating point:
      frndnum=double(ivn)*basein;
      return(1.-frndnum);
}

//==========================================================
//  initialize from a random seed iran, 0 <= iran < 259200, 
//  with the help of a "quick and dirty generator
//

void rini(int iran){


      int ir,is,irs,ir1,ikeep,ithrow,ivl1,ivl;
	  int im,ia,ic,jran,i,ifac;
      long ibase;
	  double basein;

		ir=24;
		is=10;
		irs=ir-is;
		ir1=ir-1;
		ikeep=ir;
		ithrow=24;
		ivl1=ikeep+ithrow;
		ivl=ivl1+ir1;
		ibase=long(pow(2,24));
		basein=1.0/ibase;
		srand( iran );	
		printf("Random seed: %d \n",iran);
		//iran=int(iran*quick_and_dirty_rand()+quick_and_dirty_rand());
		iran = int(quick_and_dirty_rand()*259200);
	
//  parameters for auxiliary generator used for initialization:
		
		
	im=259200;ia=7141;ic=54773;

	if((iran<0) || (iran>=im)) iran=im-1;

	jran=iran;
//  warm up the auxiliary generator a little
	for (i=1;i<=10;i++)
// cycle of auxiliary generator:
        jran=(jran*ia+ic)%im;
	
      ifac=(ibase-1)/im;
//  initialize 0 <= iv < ibase:
		for(i=0;i<=ir1;i++){
			jran=(jran*ia+ic)%im;
			iv[i]=ifac*jran;
		}

      icall=0;
}
