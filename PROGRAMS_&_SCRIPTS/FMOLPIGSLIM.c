/* FMOLPIGSLIM.c (14/02/2022) Includes Yang1-2 and Fhom */

/* ***************************************************** */

#include "libhdr"
#include "ranlib.h"
#define NN 2001  /* max number of NIND */
#define MM 8001  /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30

int marker, CHROM, NCRO, NLOCI, REPS, rep;
int lastinrecpoissontable;
int x, i, j, k, l, c, neutral;
int RM[31], initialgen[MM][31], p1, p2, ran_i;
int NINDNP, gmNP[NN][MM][2], gm[NN][MM][2], SEGRE[MM][31];
int initialgenNP[MM][31], chromNP[MM][31], chrom[MM][31];
unsigned long long int posNP[MM][31], pos[MM][31];
int father[NN], mother[NN], qtl_l, TOTLOCI;

int code[NN], dad[NN], mum[NN], coh[NN];
int numSNP, numSNPF, numSNPgen0, numQTLgen0, numSNPseg, numQTLseg, genBP;
int nind, NSEGLOCNP, current;
int locus[NN][MM][2], INITIALFREQ05;

double sNP[MM][31], atNP[MM][31], hsNP[MM][31], hatNP[MM][31], qNP[MM][31];
double pm_s[NN], pm_a[NN], L;
double AA, Aa, aa, qBP[MM][31], q[25][MM][31];
double s[MM][31], at[MM][31], hs[MM][31], hat[MM][31], LEQ_NP;
double mat[NN][NN];
double recpoissontable[maxmpt];
double Wc[NN], AVE_Wc;
double sum_G, sum_F, sum_F2, sum_FG, var_F;

double fpedBP[NN], fibdBP[NN], fvr1BP[NN], fvr2BP[NN], fyang1BP[NN], fyang2BP[NN], fLH1BP[NN], fLH2BP[NN], fhomBP[NN], num[25];
double EHomBP, EHetBP, HomBP[NN], sum_VR1BP[NN], sum_VR2BP[NN], sum_YANG1BP[NN], sum_YANG2BP[NN], sum_LH2BP[NN];

struct acc FpedBP[25], FibdBP[25], Fvr1BP[25], Fvr2BP[25], Fyang1BP[25], Fyang2BP[25], FLH1BP[25], FLH2BP[25], FhomBP[25], AVE_q[25];
struct acc AVE_FpedBP[25], AVE_FibdBP[25], AVE_Fvr1BP[25], AVE_Fvr2BP[25], AVE_Fyang1BP[25], AVE_Fyang2BP[25], AVE_FLH1BP[25], AVE_FLH2BP[25], AVE_FhomBP[25];
struct acc AVE_VFpedBP[25], AVE_VFibdBP[25], AVE_VFvr1BP[25], AVE_VFvr2BP[25], AVE_VFyang1BP[25], AVE_VFyang2BP[25], AVE_VFLH1BP[25], AVE_VFLH2BP[25], AVE_VFhomBP[25];

struct acc IDFpedBP, IDFibdBP, IDFvr1BP, IDFvr2BP, IDFyang1BP, IDFyang2BP, IDFLH1BP, IDFLH2BP, IDFhomBP;
struct acc AVE_IDFpedBP, AVE_IDFibdBP, AVE_IDFvr1BP, AVE_IDFvr2BP, AVE_IDFyang1BP, AVE_IDFyang2BP, AVE_IDFLH1BP, AVE_IDFLH2BP, AVE_IDFhomBP;

struct acc AVE_FpedFibdBP, AVE_FpedFvr1BP, AVE_FpedFvr2BP, AVE_FpedFyang1BP, AVE_FpedFyang2BP, AVE_FpedFLH1BP, AVE_FpedFLH2BP, AVE_FpedFhomBP;
struct acc AVE_FibdFvr1BP, AVE_FibdFvr2BP, AVE_FibdFyang1BP, AVE_FibdFyang2BP, AVE_FibdFLH1BP, AVE_FibdFLH2BP, AVE_FibdFhomBP;
struct acc AVE_Fvr1Fvr2BP, AVE_Fvr1Fyang1BP, AVE_Fvr1Fyang2BP, AVE_Fvr1FLH1BP, AVE_Fvr1FLH2BP, AVE_Fvr1FhomBP;
struct acc AVE_Fvr2Fyang1BP, AVE_Fvr2Fyang2BP, AVE_Fvr2FLH1BP, AVE_Fvr2FLH2BP, AVE_Fvr2FhomBP;
struct acc AVE_Fyang1Fyang2BP, AVE_Fyang1FLH1BP, AVE_Fyang1FLH2BP, AVE_Fyang1FhomBP;
struct acc AVE_Fyang2FLH1BP, AVE_Fyang2FLH2BP, AVE_Fyang2FhomBP;
struct acc AVE_FLH1FLH2BP, AVE_FLH1FhomBP;
struct acc AVE_FLH2FhomBP;

struct covacc FpedFibdBP, FpedFvr1BP, FpedFvr2BP, FpedFyang1BP, FpedFyang2BP, FpedFLH1BP, FpedFLH2BP, FpedFhomBP;
struct covacc FibdFvr1BP, FibdFvr2BP, FibdFyang1BP, FibdFyang2BP, FibdFLH1BP, FibdFLH2BP, FibdFhomBP;
struct covacc Fvr1Fvr2BP, Fvr1Fyang1BP, Fvr1Fyang2BP, Fvr1FLH1BP, Fvr1FLH2BP, Fvr1FhomBP;
struct covacc Fvr2Fyang1BP, Fvr2Fyang2BP, Fvr2FLH1BP, Fvr2FLH2BP, Fvr2FhomBP;
struct covacc Fyang1Fyang2BP, Fyang1FLH1BP, Fyang1FLH2BP, Fyang1FhomBP;
struct covacc Fyang2FLH1BP, Fyang2FLH2BP, Fyang2FhomBP;
struct covacc FLH1FLH2BP, FLH1FhomBP;
struct covacc FLH2FhomBP;

double NeFpedBP, NeFibdBP, NeFvr1BP, NeFvr2BP, NeFyang1BP, NeFyang2BP, NeFLH1BP, NeFLH2BP, NeFhomBP;
struct acc AVE_NeFpedBP, AVE_NeFibdBP, AVE_NeFvr1BP, AVE_NeFvr2BP, AVE_NeFyang1BP, AVE_NeFyang2BP, AVE_NeFLH1BP, AVE_NeFLH2BP, AVE_NeFhomBP;

void recombination_masks();

FILE *fptr, *fdat, *fpop, *flistIniSS, *flistSS, *flistSQ, *fped, *fmap, *fphe, *fmolout, *fgpig, *fcoh, *foutMF, *foutVF, *foutINDF_23, *foutINDF_BP_23, *foutID, *foutNe, *foutR, *foutQBP, *foutQ23;

/* ***************************************************** */

main()
{
	fptr = fopen ("dfilename.dat","w");
	flistIniSS = fopen ("list_initial_segregating_snps","w");
	flistSS = fopen ("list_segregating_snps","w");
	flistSQ = fopen ("list_segregating_qtls","w");
	foutINDF_23 = fopen ("list_individual_Fs_23","w");
	foutINDF_BP_23 = fopen ("list_individual_Fs_BP_23","w");

	foutMF = fopen ("outfileF.dat","w");
	foutVF = fopen ("outfileVF.dat","w");
	foutID = fopen ("outfileID.dat","w");
	foutNe = fopen ("outfileNe.dat","w");
	foutR = fopen ("outfileR.dat","w");
	foutQBP = fopen ("outfileQBP.dat","w");
	foutQ23 = fopen ("outfileQ23.dat","w");

	fped = fopen ("data.ped","w");
	fmap = fopen ("data.map","w");
	fphe = fopen ("data.phe","w");
	fmolout = fopen ("Fmolout.dat","w");

	getinputs();
	recombination_masks(RM);
	code_pigs();
	natural_population();

	for(rep=1; rep<=REPS; rep++)
	{
		sample();
		mating();
		frequency_genes();
		if (rep == 1)	coancestry_matrix();
		estimates_of_F();
		inbreeding_depression();
		correlations();
		estimates_Ne();
		if (rep == 1)	plink_files();
		settozero();
	}
	printout();
	writeseed();
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();
	getintandskip("NINDNP (max 2000):",&NINDNP,1,2000);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	getintandskip("CHROM (min 1, max 50):",&CHROM,1,50);
	NLOCI = 30;
	getintandskip("Neutral (1) or not (0) :",&neutral,0,1);
	getintandskip("genBP :",&genBP,0,23);
	getintandskip("current :",&current,0,1);
	getintandskip("INITIALFREQ05 (yes 1 or not 0) :",&INITIALFREQ05,0,1);
	getintandskip("REPS :",&REPS,0,infinity);
}

/* **************************************************** */

void recombination_masks (RM)
int RM[];
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);

    /* POISSON RECOMBINATION NUMBER */

    if ( (exp(-(double)L) != 0.0)&&((double)L < normalthreshold) )
    generatepoissontable((double)L, &lastinrecpoissontable, recpoissontable, maxmpt-1);
}

/* ***************************************************** */

code_pigs ()
{
	/* ***** take codes, father and mother of 1206 individuals ***** */

	fgpig=fopen("GenealogiaGuadyerbas","r");

	for (i=0; i<1206; i++)
	{
		fscanf(fgpig,"%d", &x);
		code[i] = x;
		fscanf(fgpig,"%d", &x);
		dad[i] = x;
		fscanf(fgpig,"%d", &x);
		mum[i] = x;
	//	if (tracelevel!=0)	if ((i<50)||(i>1180)) fprintf(fptr,"i=%d code=%d dad=%d mum=%d\n", i, code[i], dad[i], mum[i]);
	//	if (tracelevel!=0)	fprintf(fptr,"i=%d code=%d dad=%d mum=%d\n", i, code[i], dad[i], mum[i]);
	}

	fclose(fgpig);

	/* ***** take codes of 1206 genotyped individuals and cohort ***** */

	fcoh=fopen("Animales_genotipados_1206","r");

	for (j=0; j<1206; j++)
	{
		fscanf(fcoh,"%d", &x);
		fscanf(fcoh,"%d", &x);
		coh[j] = x;
//		if (tracelevel!=0)	if ((j<50)||(j>170)) fprintf(fptr,"j=%d cohort=%d\n", j, coh[j]);
//		if (tracelevel!=0)	fprintf(fptr,"j=%d cohort=%d\n", j, coh[j]);
	}

	for (j=0; j<1206; j++)
	if (coh[j] >= genBP)
	{
		nind ++;
	}

	/* ***** Father and mother of each individual ***** */

	for (i=31; i<1206; i++)
	{
		for (j=0; j<i; j++)
		{
			if (dad[i] == code[j])	father[i] = j;
			if (mum[i] == code[j])	mother[i] = j;
		}
//		if (tracelevel!=0)	if ((i<40)||(i>1170)) fprintf(fptr,"i=%d %d mother=%d\n", i, father[i], mother[i]);
		if (tracelevel!=0)	if (coh[i] <= 2) fprintf(fptr,"i=%d father=%d mother=%d\n", i, father[i], mother[i]);
	}

	for (i=0; i<1206; i++)
	{
//		if (tracelevel!=0)	if ((i<920)||(i>1170)) fprintf(fptr,"i=%d cohort=%d\n", i, coh[i]);
//		if (tracelevel!=0)	fprintf(fptr,"i=%d cohort=%d\n", i, coh[i]);
	}

	fclose(fcoh);
}

/* ***************************************************** */

natural_population ()
{
	int ds, g0, g1;
	unsigned long long int dpos;
	
	double dps, da, dh, dq;

	/* ***** take effects of genes ***** */

	fdat=fopen("list_allsnps","r");

	fscanf(fdat,"%d", &x);
	numSNP = x;

	NCRO = numSNP/NLOCI;
	TOTLOCI = NCRO * NLOCI;

	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%d%llu%lf%lf%lf%lf", &ds, &dpos, &dps, &da, &dh, &dq);
		chromNP[k][l] = ds;
		//fprintf(fptr,"chromNP[%d][%d]=%d\n", k, l, chromNP[k][l]);
		if (dps < -1.0) dps=(-1.0);
		if (da == -99.0) da=0.0;
 		if (da > 0.0) da=(-da);
		sNP[k][l] = dps;
		atNP[k][l] = da;
		posNP[k][l] = dpos;
		hsNP[k][l] = dh;
		hatNP[k][l] = dh;
//		if((k==19)&&(l==23)) fprintf(fptr,"k=%d l=%d posNP=%llu chromNP=%d\n", k, l, posNP[k][l], chromNP[k][l]);

		initialgenNP[k][l] = 0;
	}

	fclose(fdat);

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("dataBP.ped","r");

	for (i=0; i<NINDNP; i++)
	{
		lookfortext("IND");

		for (j=1; j<=5; j++)	fscanf(fpop,"%d", &x);

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			fscanf(fpop,"%d%d", &g0, &g1);

			if (g0 == 2)	gmNP[i][k][0]=(gmNP[i][k][0] | RM[l]);
			if (g1 == 2)	gmNP[i][k][1]=(gmNP[i][k][1] | RM[l]);
		}
	}

	fclose(fpop);

//	COMMENT for out
//	printf("\n");
//	for (i=0; i<NINDNP; i++)
//	{
//		if ((gmNP[i][0][0] & RM[5])==RM[5])	printf("1 ");
//		else						printf("0 ");
//		if ((gmNP[i][0][1] & RM[5])==RM[5])	printf("1 ");
//		else						printf("0 ");
//	}

	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))	aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		if (qNP[k][l] > 0.0)	NSEGLOCNP ++;
	} 

//	COMMENT for out
//	printf("\n LEQ_NP = %f", LEQ_NP);
//	for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)
//	printf("\n k=%d l=%d   sNP=%f  hNP=%f  qNP=%f  inigen=%d", k, l, sNP[k][l], hsNP[k][l], qNP[k][l], initialgenNP[k][l]);
}

/* ***************************************************** */

sample ()
{
	int g, chosen[NN];
	
	/* ***** sample the first 31 individuals from the Base Population ***** */

	for (i=0; i<NINDNP; i++)		chosen[i] = 0;

	for (i=0; i<31; i++)
	{
		do
		{ ran_i = (int)(uniform() * NINDNP); }
		while (chosen[ran_i] == 1);

		chosen[ran_i] = 1;
//		fprintf(fptr,"%d ran_i=%d\n", i, ran_i);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;

//			if (tracelevel!=0)	fprintf (fptr," i = %d  k = %d  gmNP0 = %d   gmNP1 = %d\n", i, k, gmNP[i][k][0], gmNP[i][k][1]);
		}
	}
	for (i=0; i<31; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];

//		if (tracelevel!=0)	fprintf (fptr," %d    gm0 = %d   gm1 = %d\n", i, gm[i][k][0], gm[i][k][1]);
	}

/*	if (tracelevel!=0)
	{
		for (k=0; k<NCRO; k++)
		{
			for (l=0; l<NLOCI; l++)
			{
				if(gm[0][k][l]==RM[l])  	fprintf (fptr,"1 ");
				else			   		fprintf (fptr,"0 ");
			}
			fprintf (fptr,"\n");
		}
	}
*/
	/* ***** take effects of genes from Base Population ***** */

	if (rep == 1)
	for (k=0; k<NCRO; k++) 
	for (l=0; l<NLOCI; l++)
	{
		chrom[k][l] = chromNP[k][l];
		//fprintf(fptr,"chrom[%d][%d]=%d\n", k, l, chrom[k][l]);	
		pos[k][l] = posNP[k][l];
		s[k][l] = sNP[k][l];
		at[k][l] = atNP[k][l];
		hs[k][l] = hsNP[k][l];
		hat[k][l] = hatNP[k][l];
//		if (tracelevel!=0)	if((k==19)&&(l==23)) fprintf(fptr,"k=%d l=%d pos=%llu chrom=%d\n", k, l, pos[k][l], chrom[k][l]);
		
		initialgen[k][l] = initialgenNP[k][l];

//		if (tracelevel!=0)    if ((k<5)&&(l<5)) fprintf(fptr,"\n k=%d l=%d s=%f hs=%f", k, l, s[k][l], hs[k][l]);
	}
}

/* ***************************************************** */

mating ()
{
	int EE[MM], FF[MM], last;
	double rnd;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;

	for (i=0; i<1206; i++)
	{
		// MULTILLELIC GENES

		for (k=0; k<NCRO; k++)
		{
			if (coh[i] <= genBP)
			{
				locus[i][k][0] = i;
				locus[i][k][1] = i+10000;
			}
			else
			{
				locus[i][k][0] = 0;
				locus[i][k][1] = 0;
			}
		}
	}

	for (i=31; i<1206; i++)
	{
		generahijo: /* ***** */;

		pm_s[i] = 1.0;
		pm_a[i] = 0.0;

		for (j=0; j<i; j++)
		{
			if (dad[i] == code[j])	p1 = j;
			if (mum[i] == code[j])	p2 = j;
		}
//		if (tracelevel!=0)   fprintf (fptr,"i=%d p1=%d p2=%d\n", i, p1, p2);

		if(L==99.0)
		{	    /* ******************* Free recombination ******************* */
			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	gm[i][k][0]=((EE[k]&gm[p1][k][0])|(FF[k]&gm[p1][k][1]));
			}
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], gm[p1][0][0], gm[p1][0][1], gm[p1][1][0], gm[p1][1][1], gm[p1][2][0], gm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	gm[i][k][1]=((EE[k]&gm[p2][k][0])|(FF[k]&gm[p2][k][1]));
			}
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, EE[0], EE[1], EE[2], gm[p1][0][0], gm[p1][0][1], gm[p1][1][0], gm[p1][1][1], gm[p1][2][0], gm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);

			// MULTIALLELIC GENES
			for (k=0; k<NCRO; k++)
			if (coh[i] > genBP)
			{
				if (uniform() < 0.5)	locus[i][k][0] = locus[p1][k][0];
				else					locus[i][k][0] = locus[p1][k][1];

				if (uniform() < 0.5)	locus[i][k][1] = locus[p2][k][0];
				else					locus[i][k][1] = locus[p2][k][1];
			}
		}
		else
		{	    /* ************** Restricted recombination ***************** */

			/* ****** Chromosome from father ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// SEGREGATION CHROMOSOMES
 
//			if ((tracelevel!=0)&&(i==31))	fprintf (fptr,"SEGREGATION OF CHROMOSOMES\n");
			last = 1;
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (chrom[k][l] != last)
				{
					last = chrom[k][l]; 
					//fprintf(fptr,"k=%d l=%d last=%d\n", k, l, last);
				 	ncrorec[k] = 1;
					pointrec[k][l] = 1;
				}
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			//if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs); 

			//if ((tracelevel!=0)&&(i==31))	fprintf (fptr,"RECOMBINATIONS\n");
			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
				//if ((tracelevel!=0)&&(i==31))	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gm[i][k][0]=((EE[k]&gm[p1][k][0])|(FF[k]&gm[p1][k][1]));
			}

			// MULTIALLELIC GENES
			if (coh[i] > genBP)
			for (k=0; k<NCRO; k++)
			{
				if ((EE[k] & RM[0]) == RM[0])	locus[i][k][0] = locus[p1][k][0];
				else						locus[i][k][0] = locus[p1][k][1];
			}

/*			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l]) fprintf (fptr,"1");
						else			   		fprintf (fptr,"0");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l]) fprintf (fptr,"1");
						else			   		fprintf (fptr,"0");
					}
					fprintf (fptr,"\n");
				}
			}
*/
			/* ****** Chromosome from mother ****** */

			for (k=0; k<NCRO; k++)
			{
				ncrorec[k] = 0;
				for (l=0; l<NLOCI; l++)  pointrec[k][l] = 0;
			}

			// SEGREGATION CHROMOSOMES
 
			last = 1;
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (chrom[k][l] != last)
				{
					last = chrom[k][l]; 
				 	ncrorec[k] = 1;
					pointrec[k][l] = 1;
				}
			}

			// CROSSINGOVERS

			numberrecs = recombinationnumber();
			//if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs); 

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
				//if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

			marker = 1;

			for (k=0; k<NCRO; k++)
			{
				EE[k]=0;
				if (ncrorec[k] == 0)
				{
					if (marker==(-1))
					{
						EE[k] = ~EE[k];
					}
				}
				else
				{
					for (l=0; l<NLOCI; l++)
			      	{
						if (pointrec[k][l] == 0)
						{
							if (marker==(-1))  EE[k] = EE[k] | RM[l];
						}
						else
						{
							if (marker==1)
							{
								EE[k] = EE[k] | RM[l];
								marker = marker * (-1);
							}
							else
							{
								marker = marker * (-1);
							}
						}
					}
				}
			}

			rnd = uniform();
			for (k=0; k<NCRO; k++)
			{
				if (rnd < 0.5)	EE[k] = ~EE[k];
				FF[k] = ~EE[k];
				gm[i][k][1]=((EE[k]&gm[p2][k][0])|(FF[k]&gm[p2][k][1]));
			}

			// MULTIALLELIC GENES
			for (k=0; k<NCRO; k++)
			if (coh[i] > genBP)
			{
				if ((EE[k] & RM[0]) == RM[0])	locus[i][k][1] = locus[p2][k][0];
				else						locus[i][k][1] = locus[p2][k][1];
			}

//			if (tracelevel!=0)
//			for (k=0; k<NCRO; k++)
//				if(k==0)  fprintf(fptr,"locus c=%d p1=%d p2=%d i=%d k=%d %d %d\n", coh[i], p1, p2, i, k, locus[i][k][0], locus[i][k][1]);

/*			if (tracelevel!=0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l]) fprintf (fptr,"1");
						else			   		fprintf (fptr,"0");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l]) fprintf (fptr,"1");
						else			   		fprintf (fptr,"0");
					}
					fprintf (fptr,"\n");
				}
			}
*/
		}

		/* *****Genotypic value of offspring***** */

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_a[i] += at[k][l];
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 		/* AA */;
			else	pm_a[i] += (at[k][l]*hat[k][l]);
		}

		/* *****Fitness value of offspring***** */

    		if (neutral == 0)
		{
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= (1.0 + s[k][l]);
				else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		else	pm_s[i] = 1.0;
		//if (tracelevel!=0)	fprintf (fptr," %d   pm_s = %f\n", i, pm_s[i]);

		if (uniform() > pm_s[i]) goto generahijo;
	}
}

/* ***************************************************** */

int recombinationnumber ()
{
	int r;
	if ((L < normalthreshold) && (exp(-L) != 0.0) )
	{
		r = poisson(lastinrecpoissontable, recpoissontable);
	}
	else r = (int)normal(L, sqrt(L));
	return(r);
}

/* ***************************************************** */

frequency_genes ()
{
	for (i=0; i<1206; i++)	Wc[i] = 1.0;

//	FITNESS
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		for (i=0; i<1206; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	Wc[i] *= (1.0 + s[k][l]);
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 		/* AA */;
			else	Wc[i] *= (1.0 + (s[k][l]*hs[k][l]));
		}
	}
	for (i=0; i<1206; i++)
	{
		AVE_Wc += Wc[i]/1206.0;
//		if (tracelevel!=0)	fprintf (fptr,"Wc[%d] = %f\n", i, Wc[i]);
	}
//	if (tracelevel!=0)	fprintf (fptr,"AVE_Wc = %f\n", AVE_Wc);

//	FREQUENCIES
	if (rep == 1)
	{
		for (i=0; i<1206; i++)
		for (c=0; c<=23; c++)	if (coh[i] == c)	num[c] ++;


//		if (tracelevel!=0)	for (i=0; i<100; i++)	fprintf (fptr,"\n i=%d coh=%d", i, coh[i]);
		if (tracelevel!=0)	for (c=0; c<=23; c++)	fprintf (fptr,"\n c=%d num=%f", c, num[c]);
		if (tracelevel!=0)	fprintf (fptr,"\n");

		numSNPgen0 = 0;
		numQTLgen0 = 0;
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		for (c=0; c<=23; c++)
		{
			AA=0.0; Aa=0.0; aa=0.0;

			for (i=0; i<1206; i++)
			{
				if (coh[i] == c)
				{
					if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
			    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     		else	Aa+=1.0;
				}
			}
			q[c][k][l] = (aa/num[c])+(Aa/(2.0*num[c]));
			accum(&AVE_q[c], q[c][k][l]);
//			if (tracelevel!=0)	if((k==19)&&(l==23)&&(c>=genBP)) fprintf(fptr,"k=%d l=%d chrom=%d q=%f\n", k, l, chrom[k][l], q[c][k][l]);
			if ((rep == 1)&&(c==genBP)&&(q[genBP][k][l]!=0.0)&&(q[genBP][k][l]!=1.0))	fprintf(foutQBP, "%d %f\n", chrom[k][l], q[genBP][k][l]);
			if ((rep == 1)&&(c==23)&&(q[genBP][k][l]!=0.0)&&(q[genBP][k][l]!=1.0))	fprintf(foutQ23, "%d %f\n", chrom[k][l], q[23][k][l]);
		}

		if ((q[0][k][l]!=0.0)&&(rep == 1))
		{
			numSNPgen0 ++;
			fprintf(flistIniSS,"%d %llu %f %f %f %f %f\n", chrom[k][l], pos[k][l], s[k][l], at[k][l], hs[k][l], hat[k][l], q[0][k][l]);
			if (at[k][l] != 0.0) numQTLgen0 ++;
		}
	}

//	for (k=0; k<NCRO; k++)
//	for (l=0; l<NLOCI; l++)
//	for (c=0; c<=23; c++)
//	if ((k==19)&&(l==23)&&(c>=genBP))
//	if (tracelevel!=0)	fprintf (fptr,"\n k=%d l=%d coh=%d num=%f chrom=%d q=%f", k, l, c, num[c], chrom[k][l], q[c][k][l]);

	if (rep == 1)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)	SEGRE[k][l] = 0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		for (c=genBP; c<=23; c++)
		{
			if (q[c][k][l] != 0.0)	SEGRE[k][l] = 1;
//			if (tracelevel!=0)	if((k==19)&&(l==23)&&(c>=genBP)) fprintf(fptr,"\nk=%d l=%d chrom=%d q=%f SEGRE=%d\n", k, l, chrom[k][l], q[c][k][l], SEGRE[k][l]);
		}

		numSNPseg = 0;
		numQTLseg = 0;

		//SEGREGATING SNPs
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (SEGRE[k][l] == 1)
		{
			numSNPseg ++;
			fprintf(flistSS,"%d %llu %f %f %f %f\n", chrom[k][l], pos[k][l], s[k][l], at[k][l], hs[k][l], hat[k][l]);
		}

		//SEGREGATING QTLs
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((s[k][l] != 0.0)&&(SEGRE[k][l] == 1))
		{
			numQTLseg ++;
			fprintf(flistSQ,"%d %llu %f %f %f %f\n", chrom[k][l], pos[k][l], s[k][l], at[k][l], hs[k][l], hat[k][l]);
		}
	}
}

/* ***************************************************** */

coancestry_matrix()
{
	for (i=0; i<1206; i++)
	for (j=i; j<1206; j++)
	{
		if ((coh[i] <= genBP) && (coh[j] <= genBP))
		{
			if (i == j)	mat[i][i] = 0.5;
			else			mat[i][j] = 0.0;
		}
		else
		{
			if (i == j)	mat[i][i] = 0.5 * (1.0 + mat[father[i]][mother[i]]);
			else			
			{
				mat[i][j] = 0.5 * (mat[i][father[j]] + mat[i][mother[j]]);
				mat[j][i] = mat[i][j];
			}
		}
	}

	if (tracelevel!=0)
	{
		fprintf(fptr, "\n\ni=93  coh=%d father_93=%d mother_93=%d dad_93=%d mum_93=%d\n\n", coh[93], father[93], mother[93], dad[93], mum[93]);
		fprintf(fptr, "mat_f93=%f\n\n", mat[father[93]][mother[93]]);

		fprintf(fptr, "\n\ni=69  coh=%d father_69=%d mother_69=%d dad_69=%d mum_69=%d\n\n", coh[69], father[69], mother[69], dad[69], mum[69]);
		fprintf(fptr, "mat_f69=%f\n\n", mat[father[69]][mother[69]]);

		fprintf(fptr, "\n\ni=20  coh=%d father_20=%d mother_20=%d dad_20=%d mum_20=%d\n\n", coh[20], father[20], mother[20], dad[20], mum[20]);
		fprintf(fptr, "mat_f20=%f\n\n", mat[father[20]][mother[20]]);
	}
	if (tracelevel!=0)    
	for (i=0; i<1206; i++)
	if (coh[i] <= 2)
	{
		fprintf(fptr, "i=%d   coh=%d   fat=%d   mot=%d   mat=%f   Fp=%f\n", i, coh[i], father[i], mother[i], mat[father[20]][mother[20]], 2.0*mat[i][i]-1.0);
	}
}

/* ***************************************************** */

estimates_of_F ()
{
	numSNPF = 0;
	EHomBP = 0.0;
	EHetBP = 0.0;

	if (current == 0)
	{
		// INITIAL GENERATION genBP

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)	qBP[k][l] = q[genBP][k][l]; 
	}
	else
	{
		// FREQUENCIES OF GENERATION 23

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)	qBP[k][l] = q[23][k][l]; 
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qBP[k][l] != 0.0)&&(qBP[k][l] != 1.0))
	{
		numSNPF ++;
		if (INITIALFREQ05 == 1)
		{
			EHomBP += 0.5;
			EHetBP += 0.5;
		}
		else
		{
			EHomBP += (1.0 - 2.0*qBP[k][l]*(1.0-qBP[k][l]));
			EHetBP += 2.0*qBP[k][l]*(1.0-qBP[k][l]);
		}
	}

	for (c=0; c<=23; c++)
	{
		for (i=0; i<1206; i++)
		if (coh[i] == c)
		{
			HomBP[i] = 0.0;
			sum_VR1BP[i] = 0.0;
			sum_VR2BP[i] = 0.0; 
			sum_YANG1BP[i] = 0.0; 
			sum_YANG2BP[i] = 0.0; 
			sum_LH2BP[i] = 0.0; 

			// MULTILLELIC GENES

			fibdBP[i] = 0.0;
			for (k=0; k<NCRO; k++)
			{
				if (locus[i][k][0] == locus[i][k][1])
				fibdBP[i] ++;
			}
			fibdBP[i] = fibdBP[i] / NCRO;

			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			if ((qBP[k][l] != 0.0)&&(qBP[k][l] != 1.0))
			{
//				if (tracelevel!=0) if((k==0)&&(l<5))   fprintf(fptr,"\n q[%d][%d]=%f   EHomBP=%f", k, l, qBP[k][l], EHomBP);

				if (INITIALFREQ05 == 1)
				{
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
					{
						HomBP[i] ++;
						sum_VR1BP[i] += (2.0-(2.0*0.5))*(2.0-(2.0*0.5));
						sum_VR2BP[i] += ( ( (2.0-2.0*0.5)*(2.0-2.0*0.5) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
						sum_YANG1BP[i] += ( (4.0)-((1.0+2.0*0.5)*2.0)+(2.0*0.5*0.5) ); 
						sum_YANG2BP[i] += ( (4.0)-((1.0+2.0*0.5)*2.0)+(2.0*0.5*0.5) ) / (2.0*0.5*(1.0-0.5)); 
						}
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
					{
						HomBP[i] ++;
						sum_VR1BP[i] += (0.0-(2.0*0.5))*(0.0-(2.0*0.5));
						sum_VR2BP[i] += ( ( (0.0-(2.0*0.5))*(0.0-(2.0*0.5)) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
						sum_YANG1BP[i] += (2.0*0.5*0.5); 
						sum_YANG2BP[i] += (2.0*0.5*0.5) / (2.0*0.5*(1.0-0.5)); 
					}
					else
					{
						sum_VR1BP[i] += (1.0-(2.0*0.5))*(1.0-(2.0*0.5));
						sum_VR2BP[i] += ( ( (1.0-(2.0*0.5))*(1.0-(2.0*0.5)) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
						sum_YANG1BP[i] += ( (1.0)-((1.0+2.0*0.5)*1.0)+(2.0*0.5*0.5) ); 
						sum_YANG2BP[i] += ( (1.0)-((1.0+2.0*0.5)*1.0)+(2.0*0.5*0.5) ) / (2.0*0.5*(1.0-0.5)); 
						sum_LH2BP[i] += ( 1.0 / (2.0*0.5*(1.0-0.5)) ); 
					}
				}
				else
				{
			    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
					{
						HomBP[i] ++;
						sum_VR1BP[i] += (2.0-(2.0*qBP[k][l]))*(2.0-(2.0*qBP[k][l]));
						sum_VR2BP[i] += ( ( (2.0-2.0*qBP[k][l])*(2.0-2.0*qBP[k][l]) ) / (2.0*qBP[k][l]*(1.0-qBP[k][l])) ) - 1.0; 
						sum_YANG1BP[i] += ( (4.0)-((1.0+2.0*qBP[k][l])*2.0)+(2.0*qBP[k][l]*qBP[k][l]) ); 
						sum_YANG2BP[i] += ( (4.0)-((1.0+2.0*qBP[k][l])*2.0)+(2.0*qBP[k][l]*qBP[k][l]) ) / (2.0*qBP[k][l]*(1.0-qBP[k][l])); 
					}
					else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
					{
						HomBP[i] ++;
						sum_VR1BP[i] += (0.0-(2.0*qBP[k][l]))*(0.0-(2.0*qBP[k][l]));
						sum_VR2BP[i] += ( ( (0.0-(2.0*qBP[k][l]))*(0.0-(2.0*qBP[k][l])) ) / (2.0*qBP[k][l]*(1.0-qBP[k][l])) ) - 1.0; 
						sum_YANG1BP[i] += (2.0*qBP[k][l]*qBP[k][l]); 
						sum_YANG2BP[i] += (2.0*qBP[k][l]*qBP[k][l]) / (2.0*qBP[k][l]*(1.0-qBP[k][l])); 
					}
					else
					{
						sum_VR1BP[i] += (1.0-(2.0*qBP[k][l]))*(1.0-(2.0*qBP[k][l]));
						sum_VR2BP[i] += ( ( (1.0-(2.0*qBP[k][l]))*(1.0-(2.0*qBP[k][l])) ) / (2.0*qBP[k][l]*(1.0-qBP[k][l])) ) - 1.0; 
						sum_YANG1BP[i] += ( (1.0)-((1.0+2.0*qBP[k][l])*1.0)+(2.0*qBP[k][l]*qBP[k][l]) ); 
						sum_YANG2BP[i] += ( (1.0)-((1.0+2.0*qBP[k][l])*1.0)+(2.0*qBP[k][l]*qBP[k][l]) ) / (2.0*qBP[k][l]*(1.0-qBP[k][l])); 
						sum_LH2BP[i] += ( 1.0 / (2.0*qBP[k][l]*(1.0-qBP[k][l])) ); 
					}
				}
			}

			fpedBP[i] = 2.0*mat[i][i]-1.0;
			fvr1BP[i] = (sum_VR1BP[i]/EHetBP)-1.0;
			fvr2BP[i] = sum_VR2BP[i] / numSNPF;
			fyang1BP[i] = sum_YANG1BP[i] / EHetBP;
			fyang2BP[i] = sum_YANG2BP[i] / numSNPF;
			fLH1BP[i] = (HomBP[i]-EHomBP) / (numSNPF-EHomBP);
			fLH2BP[i] = 1.0 - ( sum_LH2BP[i] / numSNPF );
			fhomBP[i] = HomBP[i] / numSNPF;

			accum(&FpedBP[c], 2.0*mat[i][i]-1.0);
			accum(&FibdBP[c], fibdBP[i]);
			accum(&Fvr1BP[c], fvr1BP[i]);
			accum(&Fvr2BP[c], fvr2BP[i]);
			accum(&Fyang1BP[c], fyang1BP[i]);
			accum(&Fyang2BP[c], fyang2BP[i]);
			accum(&FLH1BP[c], fLH1BP[i]);
			accum(&FLH2BP[c], fLH2BP[i]);
			accum(&FhomBP[c], fhomBP[i]);

//			fprintf(fptr,"%f %f %f %f %f %f %f\n", fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
			if (coh[i]>=genBP)	fprintf(foutINDF_BP_23,"%d %d %f %f %f %f %f %f %f %f %f\n", code[i], coh[i], 2.0*mat[i][i]-1.0, fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
			if (coh[i]==23)	fprintf(foutINDF_23,"%d %d %f %f %f %f %f %f %f %f %f\n", code[i], coh[i], 2.0*mat[i][i]-1.0, fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);

/*			if (tracelevel!=0)	
			{
				fprintf(fptr,"\ni=%d  sumVR1BP=%f  sumVR2BP=%f  sumYang1BP=%f  sumYang2BP=%f  HomBP=%f  EHomBP=%f\n", i, sum_VR1BP[i], sum_VR2BP[i], sum_YANG1BP[i], sum_YANG2BP[i], HomBP[i], EHomBP);
				fprintf(fptr,"i=%d  numSNPF=%d  fvr1BP=%f  fvr2BP=%f  fyang1BP=%f  fyang2BP=%f  fLH1BP=%f  fLH2BP=%f  fhomBP=%f\n", i, numSNPF, fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
			}
			fprintf(fptr,"i=%d  pm_s=%f  numLOCI=%d  fvr1=%f  fvr2=%f  fyang1=%f  fyang2=%f  fLH1=%f  fLH2=%f  fhom=%f\n", i, pm_s[i], numSNPF, fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);
*/
		}
	}

	for (c=0; c<=23; c++)
	{
		accum(&AVE_FpedBP[c], accmean(&FpedBP[c]));
		accum(&AVE_FibdBP[c], accmean(&FibdBP[c]));
		accum(&AVE_Fvr1BP[c], accmean(&Fvr1BP[c]));
		accum(&AVE_Fvr2BP[c], accmean(&Fvr2BP[c]));
		accum(&AVE_Fyang1BP[c], accmean(&Fyang1BP[c]));
		accum(&AVE_Fyang2BP[c], accmean(&Fyang2BP[c]));
		accum(&AVE_FLH1BP[c], accmean(&FLH1BP[c]));
		accum(&AVE_FLH2BP[c], accmean(&FLH2BP[c]));
		accum(&AVE_FhomBP[c], accmean(&FhomBP[c]));

		accum(&AVE_VFpedBP[c], variance(&FpedBP[c]));
		accum(&AVE_VFibdBP[c], variance(&FibdBP[c]));
		accum(&AVE_VFvr1BP[c], variance(&Fvr1BP[c]));
		accum(&AVE_VFvr2BP[c], variance(&Fvr2BP[c]));
		accum(&AVE_VFyang1BP[c], variance(&Fyang1BP[c]));
		accum(&AVE_VFyang2BP[c], variance(&Fyang2BP[c]));
		accum(&AVE_VFLH1BP[c], variance(&FLH1BP[c]));
		accum(&AVE_VFLH2BP[c], variance(&FLH2BP[c]));
		accum(&AVE_VFhomBP[c], variance(&FhomBP[c]));

		if (c == 23)
		{
			fprintf(foutMF, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			accmean(&FpedBP[c]), accmean(&FibdBP[c]), accmean(&Fvr1BP[c]), accmean(&Fvr2BP[c]), accmean(&Fyang1BP[c]), accmean(&Fyang2BP[c]), accmean(&FLH1BP[c]), accmean(&FLH2BP[c]), accmean(&FhomBP[c]));

			fprintf(foutVF, "%f  %f  %f  %f  %f  %f  %f  %f  %f\n",
			variance(&FpedBP[c]), variance(&FibdBP[c]), variance(&Fvr1BP[c]), variance(&Fvr2BP[c]), variance(&Fyang1BP[c]), variance(&Fyang2BP[c]), variance(&FLH1BP[c]), variance(&FLH2BP[c]), variance(&FhomBP[c]));
		}
	}
}

/* ***************************************************** */

inbreeding_depression ()
{
	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fpedBP[i];
		sum_F2 += (fpedBP[i] * fpedBP[i]);
		sum_FG += (fpedBP[i] * pm_a[i]);	

		if (c == 5) fprintf(fptr,"i=%d pm=%f FpedBP=%f sum_G=%f sum_F=%f sum_F2=%f sum_FG=%f\n",
		i, pm_s[i], fpedBP[i], sum_G, sum_F, sum_F2, sum_FG); 
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);

	if (var_F > 0.0) accum (&IDFpedBP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	if (c == 5) fprintf(fptr,"\nc=%d nind=%d var_F=%f ID=%f\n", c, nind, var_F, accmean(&IDFpedBP)); 

	// *********************** FibdBP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fibdBP[i];
		sum_F2 += (fibdBP[i] * fibdBP[i]);
		sum_FG += (fibdBP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFibdBP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** Fvr1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fvr1BP[i];
		sum_F2 += (fvr1BP[i] * fvr1BP[i]);
		sum_FG += (fvr1BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFvr1BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** Fvr2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fvr2BP[i];
		sum_F2 += (fvr2BP[i] * fvr2BP[i]);
		sum_FG += (fvr2BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFvr2BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** Fyang1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fyang1BP[i];
		sum_F2 += (fyang1BP[i] * fyang1BP[i]);
		sum_FG += (fyang1BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFyang1BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** Fyang2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fyang2BP[i];
		sum_F2 += (fyang2BP[i] * fyang2BP[i]);
		sum_FG += (fyang2BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFyang2BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** FLH1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fLH1BP[i];
		sum_F2 += (fLH1BP[i] * fLH1BP[i]);
		sum_FG += (fLH1BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFLH1BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** FLH2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fLH2BP[i];
		sum_F2 += (fLH2BP[i] * fLH2BP[i]);
		sum_FG += (fLH2BP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFLH2BP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************** FhomBP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		sum_G += pm_a[i];
		sum_F += fhomBP[i];
		sum_F2 += (fhomBP[i] * fhomBP[i]);
		sum_FG += (fhomBP[i] * pm_a[i]);
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)nind);
	if (var_F > 0.0) accum (&IDFhomBP, (sum_FG - (sum_G * sum_F / (double)nind)) / var_F);

	// *********************************************************

	accum(&AVE_IDFpedBP, accmean(&IDFpedBP));
	accum(&AVE_IDFibdBP, accmean(&IDFibdBP));
	accum(&AVE_IDFvr1BP, accmean(&IDFvr1BP));
	accum(&AVE_IDFvr2BP, accmean(&IDFvr2BP));
	accum(&AVE_IDFyang1BP, accmean(&IDFyang1BP));
	accum(&AVE_IDFyang2BP, accmean(&IDFyang2BP));
	accum(&AVE_IDFLH1BP, accmean(&IDFLH1BP));
	accum(&AVE_IDFLH2BP, accmean(&IDFLH2BP));
	accum(&AVE_IDFhomBP, accmean(&IDFhomBP));

	fprintf(foutID, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&IDFpedBP), accmean(&IDFibdBP), accmean(&IDFvr1BP), accmean(&IDFvr2BP), accmean(&IDFyang1BP), accmean(&IDFyang2BP), accmean(&IDFLH1BP), accmean(&IDFLH2BP), accmean(&IDFhomBP));
}

/* ***************************************************** */

correlations ()
{
	for (i=0; i<1206; i++)
	if (coh[i] == 23)
//	if (coh[i] >= genBP)
	{
		covaccum (&FpedFibdBP, fpedBP[i], fibdBP[i]);
		covaccum (&FpedFvr1BP, fpedBP[i], fvr1BP[i]);
		covaccum (&FpedFvr2BP, fpedBP[i], fvr2BP[i]);
		covaccum (&FpedFyang1BP, fpedBP[i], fyang1BP[i]);
		covaccum (&FpedFyang2BP, fpedBP[i], fyang2BP[i]);
		covaccum (&FpedFLH1BP, fpedBP[i], fLH1BP[i]);
		covaccum (&FpedFLH2BP, fpedBP[i], fLH2BP[i]);
		covaccum (&FpedFhomBP, fpedBP[i], fhomBP[i]);
		covaccum (&FibdFvr1BP, fibdBP[i], fvr1BP[i]);
		covaccum (&FibdFvr2BP, fibdBP[i], fvr2BP[i]);
		covaccum (&FibdFyang1BP, fibdBP[i], fyang1BP[i]);
		covaccum (&FibdFyang2BP, fibdBP[i], fyang2BP[i]);
		covaccum (&FibdFLH1BP, fibdBP[i], fLH1BP[i]);
		covaccum (&FibdFLH2BP, fibdBP[i], fLH2BP[i]);
		covaccum (&FibdFhomBP, fibdBP[i], fhomBP[i]);
		covaccum (&Fvr1Fvr2BP, fvr1BP[i], fvr2BP[i]);
		covaccum (&Fvr1Fyang1BP, fvr1BP[i], fyang1BP[i]);
		covaccum (&Fvr1Fyang2BP, fvr1BP[i], fyang2BP[i]);
		covaccum (&Fvr1FLH1BP, fvr1BP[i], fLH1BP[i]);
		covaccum (&Fvr1FLH2BP, fvr1BP[i], fLH2BP[i]);
		covaccum (&Fvr1FhomBP, fvr1BP[i], fhomBP[i]);
		covaccum (&Fvr2Fyang1BP, fvr2BP[i], fyang1BP[i]);
		covaccum (&Fvr2Fyang2BP, fvr2BP[i], fyang2BP[i]);
		covaccum (&Fvr2FLH1BP, fvr2BP[i], fLH1BP[i]);
		covaccum (&Fvr2FLH2BP, fvr2BP[i], fLH2BP[i]);
		covaccum (&Fvr2FhomBP, fvr2BP[i], fhomBP[i]);
		covaccum (&Fyang1Fyang2BP, fyang1BP[i], fyang2BP[i]);
		covaccum (&Fyang1FLH1BP, fyang1BP[i], fLH1BP[i]);
		covaccum (&Fyang1FLH2BP, fyang1BP[i], fLH2BP[i]);
		covaccum (&Fyang1FhomBP, fyang1BP[i], fhomBP[i]);
		covaccum (&Fyang2FLH1BP, fyang2BP[i], fLH1BP[i]);
		covaccum (&Fyang2FLH2BP, fyang2BP[i], fLH2BP[i]);
		covaccum (&Fyang2FhomBP, fyang2BP[i], fhomBP[i]);
		covaccum (&FLH1FLH2BP, fLH1BP[i], fLH2BP[i]);
		covaccum (&FLH1FhomBP, fLH1BP[i], fhomBP[i]);
		covaccum (&FLH2FhomBP, fLH2BP[i], fhomBP[i]);
	}

	accum (&AVE_FpedFibdBP, correlation(&FpedFibdBP));
	accum (&AVE_FpedFvr1BP, correlation(&FpedFvr1BP));
	accum (&AVE_FpedFvr2BP, correlation(&FpedFvr2BP));
	accum (&AVE_FpedFyang1BP, correlation(&FpedFyang1BP));
	accum (&AVE_FpedFyang2BP, correlation(&FpedFyang2BP));
	accum (&AVE_FpedFLH1BP, correlation(&FpedFLH1BP));
	accum (&AVE_FpedFLH2BP, correlation(&FpedFLH2BP));
	accum (&AVE_FpedFhomBP, correlation(&FpedFhomBP));
	accum (&AVE_FibdFvr1BP, correlation(&FibdFvr1BP));
	accum (&AVE_FibdFvr2BP, correlation(&FibdFvr2BP));
	accum (&AVE_FibdFyang1BP, correlation(&FibdFyang1BP));
	accum (&AVE_FibdFyang2BP, correlation(&FibdFyang2BP));
	accum (&AVE_FibdFLH1BP, correlation(&FibdFLH1BP));
	accum (&AVE_FibdFLH2BP, correlation(&FibdFLH2BP));
	accum (&AVE_FibdFhomBP, correlation(&FibdFhomBP));
	accum (&AVE_Fvr1Fvr2BP, correlation(&Fvr1Fvr2BP));
	accum (&AVE_Fvr1Fyang1BP, correlation(&Fvr1Fyang1BP));
	accum (&AVE_Fvr1Fyang2BP, correlation(&Fvr1Fyang2BP));
	accum (&AVE_Fvr1FLH1BP, correlation(&Fvr1FLH1BP));
	accum (&AVE_Fvr1FLH2BP, correlation(&Fvr1FLH2BP));
	accum (&AVE_Fvr1FhomBP, correlation(&Fvr1FhomBP));
	accum (&AVE_Fvr2Fyang1BP, correlation(&Fvr2Fyang1BP));
	accum (&AVE_Fvr2Fyang2BP, correlation(&Fvr2Fyang2BP));
	accum (&AVE_Fvr2FLH1BP, correlation(&Fvr2FLH1BP));
	accum (&AVE_Fvr2FLH2BP, correlation(&Fvr2FLH2BP));
	accum (&AVE_Fvr2FhomBP, correlation(&Fvr2FhomBP));
	accum (&AVE_Fyang1Fyang2BP, correlation(&Fyang1Fyang2BP));
	accum (&AVE_Fyang1FLH1BP, correlation(&Fyang1FLH1BP));
	accum (&AVE_Fyang1FLH2BP, correlation(&Fyang1FLH2BP));
	accum (&AVE_Fyang1FhomBP, correlation(&Fyang1FhomBP));
	accum (&AVE_Fyang2FLH1BP, correlation(&Fyang2FLH1BP));
	accum (&AVE_Fyang2FLH2BP, correlation(&Fyang2FLH2BP));
	accum (&AVE_Fyang2FhomBP, correlation(&Fyang2FhomBP));
	accum (&AVE_FLH1FLH2BP, correlation(&FLH1FLH2BP));
	accum (&AVE_FLH1FhomBP, correlation(&FLH1FhomBP));
	accum (&AVE_FLH2FhomBP, correlation(&FLH2FhomBP));

	fprintf(foutR, "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
	correlation(&FpedFibdBP), correlation(&FpedFvr1BP), correlation(&FpedFvr2BP), correlation(&FpedFyang1BP), correlation(&FpedFyang2BP), correlation(&FpedFLH1BP), correlation(&FpedFLH2BP), correlation(&FpedFhomBP),
	correlation(&FibdFvr1BP), correlation(&FibdFvr2BP), correlation(&FibdFyang1BP), correlation(&FibdFyang2BP), correlation(&FibdFLH1BP), correlation(&FibdFLH2BP), correlation(&FibdFhomBP),
	correlation(&Fvr1Fvr2BP), correlation(&Fvr1Fyang1BP), correlation(&Fvr1Fyang2BP), correlation(&Fvr1FLH1BP), correlation(&Fvr1FLH2BP), correlation(&Fvr1FhomBP), correlation(&Fvr2Fyang1BP), correlation(&Fvr2Fyang2BP),
	correlation(&Fvr2FLH1BP), correlation(&Fvr2FLH2BP), correlation(&Fvr2FhomBP), correlation(&Fyang1Fyang2BP), correlation(&Fyang1FLH1BP), correlation(&Fyang1FLH2BP), correlation(&Fyang1FhomBP), correlation(&Fyang2FLH1BP), correlation(&Fyang2FLH2BP), correlation(&Fyang2FhomBP), correlation(&FLH1FLH2BP), correlation(&FLH1FhomBP), correlation(&FLH2FhomBP));
}

/* ***************************************************** */

estimates_Ne ()
{
	NeFpedBP = (23.0-genBP) / ( 2.0*(accmean(&FpedBP[23])-accmean(&FpedBP[genBP]))/(1.0-accmean(&FpedBP[genBP])) );
	NeFibdBP = (23.0-genBP) / ( 2.0*(accmean(&FibdBP[23])-accmean(&FibdBP[genBP]))/(1.0-accmean(&FibdBP[genBP])) );
	NeFvr1BP = (23.0-genBP) / ( 2.0*(accmean(&Fvr1BP[23])-accmean(&Fvr1BP[genBP]))/(1.0-accmean(&Fvr1BP[genBP])) );
	NeFvr2BP = (23.0-genBP) / ( 2.0*(accmean(&Fvr2BP[23])-accmean(&Fvr2BP[genBP]))/(1.0-accmean(&Fvr2BP[genBP])) );
	NeFyang1BP = (23.0-genBP) / ( 2.0*(accmean(&Fyang1BP[23])-accmean(&Fyang1BP[genBP]))/(1.0-accmean(&Fyang1BP[genBP])) );
	NeFyang2BP = (23.0-genBP) / ( 2.0*(accmean(&Fyang2BP[23])-accmean(&Fyang2BP[genBP]))/(1.0-accmean(&Fyang2BP[genBP])) );
	NeFLH1BP = (23.0-genBP) / ( 2.0*(accmean(&FLH1BP[23])-accmean(&FLH1BP[genBP]))/(1.0-accmean(&FLH1BP[genBP])) );
	NeFLH2BP = (23.0-genBP) / ( 2.0*(accmean(&FLH2BP[23])-accmean(&FLH2BP[genBP]))/(1.0-accmean(&FLH2BP[genBP])) );
	NeFhomBP = (23.0-genBP) / ( 2.0*(accmean(&FhomBP[23])-accmean(&FhomBP[genBP]))/(1.0-accmean(&FhomBP[genBP])) );

	accum (&AVE_NeFpedBP, NeFpedBP);
	accum (&AVE_NeFibdBP, NeFibdBP);
	accum (&AVE_NeFvr1BP, NeFvr1BP);
	accum (&AVE_NeFvr2BP, NeFvr2BP);
	accum (&AVE_NeFyang1BP, NeFyang1BP);
	accum (&AVE_NeFyang2BP, NeFyang2BP);
	accum (&AVE_NeFLH1BP, NeFLH1BP);
	accum (&AVE_NeFLH2BP, NeFLH2BP);
	accum (&AVE_NeFhomBP, NeFhomBP);

	fprintf(foutNe, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", NeFpedBP, NeFibdBP, NeFvr1BP, NeFvr2BP, NeFyang1BP, NeFyang2BP, NeFLH1BP, NeFLH2BP, NeFhomBP);
}

/* ***************************************************** */

plink_files()
{
	// data.map

	for (k=0; k<NCRO; k++)
	{
		for (l=0; l<NLOCI; l++)
		{
			if (SEGRE[k][l] == 1)
			{
				fprintf(fmap,"%d SNP%llu 0 %llu\n", chrom[k][l], pos[k][l], pos[k][l]);
//				if((k==19)&&(l==23)) fprintf(fptr,"\n%d SNP%llu 0 %llu\n", chrom[k][l], pos[k][l], pos[k][l]);
			}
		}
	}

	for (i=0; i<1206; i++)
	if (coh[i] >= genBP)
	{
		// data.phe

 		fprintf(fphe,"%d %f %f %f %f %f %f %f %f %f %f %f\n", code[i], pm_s[i], pm_a[i], fpedBP[i], fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);

		// data.ped

		fprintf(fped,"1 %d 0 0 1 -9 ", code[i]);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (SEGRE[k][l] == 1)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(fped,"T T  ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(fped,"A A  ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(fped,"T A  ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(fped,"A T  ");
		}
		fprintf(fped,"\n");
	}

	return(0);
}

/* ***************************************************** */

settozero()
{
	for (c=0; c<=23; c++)
	{
		initacc (&FpedBP[c]);
		initacc (&FibdBP[c]);
		initacc (&Fvr1BP[c]);
		initacc (&Fvr2BP[c]);
		initacc (&Fyang1BP[c]);
		initacc (&Fyang2BP[c]);
		initacc (&FLH1BP[c]);
		initacc (&FLH2BP[c]);
		initacc (&FhomBP[c]);
	}

	initacc (&IDFpedBP);
	initacc (&IDFibdBP);
	initacc (&IDFvr1BP);
	initacc (&IDFvr2BP);
	initacc (&IDFyang1BP);
	initacc (&IDFyang2BP);
	initacc (&IDFLH1BP);
	initacc (&IDFLH2BP);
	initacc (&IDFhomBP);

	initcovacc (&FpedFibdBP);
	initcovacc (&FpedFvr1BP);
	initcovacc (&FpedFvr2BP);
	initcovacc (&FpedFyang1BP);
	initcovacc (&FpedFyang2BP);
	initcovacc (&FpedFLH1BP);
	initcovacc (&FpedFLH2BP);
	initcovacc (&FpedFhomBP);
	initcovacc (&FibdFvr1BP);
	initcovacc (&FibdFvr2BP);
	initcovacc (&FibdFyang1BP);
	initcovacc (&FibdFyang2BP);
	initcovacc (&FibdFLH1BP);
	initcovacc (&FibdFLH2BP);
	initcovacc (&FibdFhomBP);
	initcovacc (&Fvr1Fvr2BP);
	initcovacc (&Fvr1Fyang1BP);
	initcovacc (&Fvr1Fyang2BP);
	initcovacc (&Fvr1FLH1BP);
	initcovacc (&Fvr1FLH2BP);
	initcovacc (&Fvr1FhomBP);
	initcovacc (&Fvr2Fyang1BP);
	initcovacc (&Fvr2Fyang2BP);
	initcovacc (&Fvr2FLH1BP);
	initcovacc (&Fvr2FLH2BP);
	initcovacc (&Fvr2FhomBP);
	initcovacc (&Fyang1Fyang2BP);
	initcovacc (&Fyang1FLH1BP);
	initcovacc (&Fyang1FLH2BP);
	initcovacc (&Fyang1FhomBP);
	initcovacc (&Fyang2FLH1BP);
	initcovacc (&Fyang2FLH2BP);
	initcovacc (&Fyang2FhomBP);
	initcovacc (&FLH1FLH2BP);
	initcovacc (&FLH1FhomBP);
	initcovacc (&FLH2FhomBP);
}

/* ***************************************************** */

printout()
{
 	fprintf(fmolout, "OUTPUT\n"); 
 	fprintf(fmolout, "NINDNP=%d NCRO=%d TOTLOCI=%d REPS=%d \nNSEGLOCNP=%d numSNPgen0=%d numQTLgen0=%d numSNPseg=%d numQTLseg=%d \nLEQ_NP=%f\n",
 		NINDNP, NCRO, TOTLOCI, REPS, NSEGLOCNP, numSNPgen0, numQTLgen0, numSNPseg, numQTLseg, LEQ_NP); 
 
 	fprintf(fmolout, "\nMEAN F\n"); 
	fprintf(fmolout, "coh   q     Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP  FhomBP\n"); 

	if (current == 0)
	{
		for (c=genBP; c<=23; c++)
		fprintf(fmolout, "%d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", c, accmean(&AVE_q[c]), accmean(&AVE_FpedBP[c]), accmean(&AVE_FibdBP[c]), accmean(&AVE_Fvr1BP[c]), accmean(&AVE_Fvr2BP[c]), accmean(&AVE_Fyang1BP[c]), accmean(&AVE_Fyang2BP[c]), accmean(&AVE_FLH1BP[c]), accmean(&AVE_FLH2BP[c]), accmean(&AVE_FhomBP[c]));
	}
	else	fprintf(fmolout, "23  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", accmean(&AVE_q[23]), accmean(&AVE_FpedBP[23]), accmean(&AVE_FibdBP[23]), accmean(&AVE_Fvr1BP[23]), accmean(&AVE_Fvr2BP[23]), accmean(&AVE_Fyang1BP[23]), accmean(&AVE_Fyang2BP[23]), accmean(&AVE_FLH1BP[23]), accmean(&AVE_FLH2BP[23]), accmean(&AVE_FhomBP[23]));

	fprintf(fmolout, "\nSD F\n"); 
	fprintf(fmolout, "coh   q     Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang1BP FLH1BP  FLH2BP  FhomBP\n"); 
	if (current == 0)
	{
		for (c=genBP; c<=23; c++)
		fprintf(fmolout, "%d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", c, sqrt(variance(&AVE_q[c])), sqrt(variance(&AVE_FpedBP[c])), sqrt(variance(&AVE_FibdBP[c])), sqrt(variance(&AVE_Fvr1BP[c])), sqrt(variance(&AVE_Fvr2BP[c])), sqrt(variance(&AVE_Fyang1BP[c])), sqrt(variance(&AVE_Fyang2BP[c])), sqrt(variance(&AVE_FLH1BP[c])), sqrt(variance(&AVE_FLH2BP[c])), sqrt(variance(&AVE_FhomBP[c])));
	}
	else fprintf(fmolout, "23  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", sqrt(variance(&AVE_q[23])), sqrt(variance(&AVE_FpedBP[23])), sqrt(variance(&AVE_FibdBP[23])), sqrt(variance(&AVE_Fvr1BP[23])), sqrt(variance(&AVE_Fvr2BP[23])), sqrt(variance(&AVE_Fyang1BP[23])), sqrt(variance(&AVE_Fyang2BP[23])), sqrt(variance(&AVE_FLH1BP[23])), sqrt(variance(&AVE_FLH2BP[23])), sqrt(variance(&AVE_FhomBP[23])));

 	fprintf(fmolout, "\nVARIANCE F\n"); 
	fprintf(fmolout, "coh   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP  FhomBP\n"); 
	if (current == 0)
	{
		for (c=genBP; c<=23; c++)
		fprintf(fmolout, "%d  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", c, accmean(&AVE_VFpedBP[c]), accmean(&AVE_VFibdBP[c]), accmean(&AVE_VFvr1BP[c]), accmean(&AVE_VFvr2BP[c]), accmean(&AVE_VFyang1BP[c]), accmean(&AVE_VFyang2BP[c]), accmean(&AVE_VFLH1BP[c]), accmean(&AVE_VFLH2BP[c]), accmean(&AVE_VFhomBP[c]));
	}
	else	fprintf(fmolout, "23  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", accmean(&AVE_VFpedBP[23]), accmean(&AVE_VFibdBP[23]), accmean(&AVE_VFvr1BP[23]), accmean(&AVE_VFvr2BP[23]), accmean(&AVE_VFyang1BP[23]), accmean(&AVE_VFyang2BP[23]), accmean(&AVE_VFLH1BP[23]), accmean(&AVE_VFLH2BP[23]), accmean(&AVE_VFhomBP[23]));

	fprintf(fmolout, "\nSD VAR F\n"); 
	fprintf(fmolout, "coh   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP  FhomBP\n"); 
	if (current == 0)
	{
		for (c=genBP; c<=23; c++)
		fprintf(fmolout, "%d  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", c, sqrt(variance(&AVE_VFpedBP[c])), sqrt(variance(&AVE_VFibdBP[c])), sqrt(variance(&AVE_VFvr1BP[c])), sqrt(variance(&AVE_VFvr2BP[c])), sqrt(variance(&AVE_VFyang1BP[c])), sqrt(variance(&AVE_VFyang2BP[c])), sqrt(variance(&AVE_VFLH1BP[c])), sqrt(variance(&AVE_VFLH2BP[c])), sqrt(variance(&AVE_VFhomBP[c])));
	}
	else fprintf(fmolout, "23  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", sqrt(variance(&AVE_VFpedBP[23])), sqrt(variance(&AVE_VFibdBP[23])), sqrt(variance(&AVE_VFvr1BP[23])), sqrt(variance(&AVE_VFvr2BP[23])), sqrt(variance(&AVE_VFyang1BP[23])), sqrt(variance(&AVE_VFyang2BP[23])), sqrt(variance(&AVE_VFLH1BP[23])), sqrt(variance(&AVE_VFLH2BP[23])), sqrt(variance(&AVE_VFhomBP[23])));

	fprintf(fmolout, "\nID\n");
	fprintf(fmolout, "nind   IDFped  IDFibd  IDFvr1  IDFvr2   IDFyang1   IDFyang2   IDFLH1   IDFLH2   IDFhom\n"); 
	fprintf(fmolout, "%3d    %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", nind, accmean(&AVE_IDFpedBP), accmean(&AVE_IDFibdBP), accmean(&AVE_IDFvr1BP), accmean(&AVE_IDFvr2BP), accmean(&AVE_IDFyang1BP), accmean(&AVE_IDFyang2BP), accmean(&AVE_IDFLH1BP), accmean(&AVE_IDFLH2BP), accmean(&AVE_IDFhomBP));

	fprintf(fmolout, "\nSD ID\n");
	fprintf(fmolout, "nind IDFped  IDFibd  IDFvr1  IDFvr2  IDFyang1 IDFyang2 IDFLH1 IDFLH2 IDFhom\n"); 
	fprintf(fmolout, "%3d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", nind, sqrt(variance(&AVE_IDFpedBP)), sqrt(variance(&AVE_IDFibdBP)), sqrt(variance(&AVE_IDFvr1BP)), sqrt(variance(&AVE_IDFvr2BP)), sqrt(variance(&AVE_IDFyang1BP)), sqrt(variance(&AVE_IDFyang2BP)), sqrt(variance(&AVE_IDFLH1BP)), sqrt(variance(&AVE_IDFLH2BP)), sqrt(variance(&AVE_IDFhomBP)));

	if (current == 0) fprintf(fmolout, "\nMEAN Correlations (at generation 23) with frequencies of cohort %d", genBP); 
	else fprintf(fmolout, "\nMEAN Correlations (at generation 23) with frequencies of cohort %d", 23); 
	fprintf(fmolout, "\n        FibdBP  Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP FLH2BP FhomBP\n"); 
	fprintf(fmolout, "Fped  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FpedFibdBP), accmean(&AVE_FpedFvr1BP), accmean(&AVE_FpedFvr2BP), accmean(&AVE_FpedFyang1BP), accmean(&AVE_FpedFyang2BP), accmean(&AVE_FpedFLH1BP), accmean(&AVE_FpedFLH2BP), accmean(&AVE_FpedFhomBP)); 
	fprintf(fmolout, "Fibd          %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FibdFvr1BP), accmean(&AVE_FibdFvr2BP), accmean(&AVE_FibdFyang1BP), accmean(&AVE_FibdFyang2BP), accmean(&AVE_FibdFLH1BP), accmean(&AVE_FibdFLH2BP), accmean(&AVE_FibdFhomBP)); 
	fprintf(fmolout, "Fvr1BP                %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr1Fvr2BP), accmean(&AVE_Fvr1Fyang1BP), accmean(&AVE_Fvr1Fyang2BP), accmean(&AVE_Fvr1FLH1BP), accmean(&AVE_Fvr1FLH2BP), accmean(&AVE_Fvr1FhomBP)); 
	fprintf(fmolout, "Fvr2BP                        %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr2Fyang1BP), accmean(&AVE_Fvr2Fyang2BP), accmean(&AVE_Fvr2FLH1BP), accmean(&AVE_Fvr2FLH2BP), accmean(&AVE_Fvr2FhomBP)); 
	fprintf(fmolout, "Fyang1BP                              %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fyang1Fyang2BP), accmean(&AVE_Fyang1FLH1BP), accmean(&AVE_Fyang1FLH2BP), accmean(&AVE_Fyang1FhomBP)); 
	fprintf(fmolout, "Fyang2BP                               	      %7.4f %7.4f %7.4f\n", accmean(&AVE_Fyang2FLH1BP), accmean(&AVE_Fyang2FLH2BP), accmean(&AVE_Fyang2FhomBP)); 
	fprintf(fmolout, "FLH1BP                                                %7.4f %7.4f\n", accmean(&AVE_FLH1FLH2BP), accmean(&AVE_FLH1FhomBP)); 
	fprintf(fmolout, "FLH2BP                                                        %7.4f\n", accmean(&AVE_FLH2FhomBP)); 

	fprintf(fmolout, "\nMEAN Ne (cohorts %d-23)", genBP); 

	fprintf(fmolout, "\nNeFped  NeFibd  NeFvr1BP NeFvr2BP  NeFyang1BP NeFyang2BP NeFLH1BP NeFLH2BP NeFhomBP\n");
	fprintf(fmolout, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_NeFpedBP), accmean(&AVE_NeFibdBP), accmean(&AVE_NeFvr1BP), accmean(&AVE_NeFvr2BP), accmean(&AVE_NeFyang1BP), accmean(&AVE_NeFyang2BP), accmean(&AVE_NeFLH1BP), accmean(&AVE_NeFLH2BP), accmean(&AVE_NeFhomBP) );

	fprintf(fmolout, "\nSD Ne");
	fprintf(fmolout, "\nNeFped  NeFibd  NeFvr1BP NeFvr2BP  NeFyang1BP NeFyang2BP NeFLH1BP NeFLH2BP NeFhomBP\n");
	fprintf(fmolout, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		sqrt(variance(&AVE_NeFpedBP)), sqrt(variance(&AVE_NeFibdBP)), sqrt(variance(&AVE_NeFvr1BP)), sqrt(variance(&AVE_NeFvr2BP)), sqrt(variance(&AVE_NeFyang1BP)), sqrt(variance(&AVE_NeFyang2BP)), sqrt(variance(&AVE_NeFLH1BP)), sqrt(variance(&AVE_NeFLH2BP)), sqrt(variance(&AVE_NeFhomBP)) );

	return(0);
}

/* ***************************************************** */

lookfortext(s)
char *s;
{
   int len, i, curchar;
   char c;

   curchar = 0;
   len = 0;

   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   do
   {
      c = getc(fpop);

      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(0);
      }
      else curchar = 0;
   }
   while (c != EOF);
}

/* ********************************************************************************************* */




