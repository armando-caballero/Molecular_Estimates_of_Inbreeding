/* EstimationMolNe3.c (including ID, Ne and artificial selection) (16/03/2022) Includes Yang1-2 y Fhom*/

/* ***************************************************** */

#include "libhdr"
#define NN 1001  /* max number of TNIND and gen */
#define MM 8001  /* max number of NCRO */
#define maxmpt 5001
#define normalthreshold 30

int ran_i, muts, lastinrecpoissontable;
int numberofrecom, marker, TNIND, NIND, NCRO, NLOCI, TOTLOCI, numLOCI[NN], numLOCIBP[NN];
int gen, genBP, generations, i, j, k, l, rep, replicates;
int RM[NN], p1, p2;
int gm[NN][MM][2], sm[NN][MM][2];
int father[MM], mother[MM];
int neutral, EC, INITIALFREQ05;
int slocus[NN][MM][2], locus[NN][MM][2];

int NINDNP, NSEGLOCNP, gmNP[NN][MM][2];
int chromNP[MM][31], chrom[MM][31];
unsigned long long int posNP[MM][31], pos[MM][31];
double sNP[MM][31], atNP[MM][31], hsNP[MM][31], hatNP[MM][31], qNP[MM][31];
double s[MM][31], at[MM][31], hs[MM][31], hat[MM][31], LEQ_NP;
double rnd, MAF;

double pm_s[NN], pm_a[NN], geno_a[NN], L, VE, va, vd, id, idgenBP, idgenLAST, ha, sd_effect;
double AA, Aa, aa, qi[MM][31], q[MM][31];
double recpoissontable[maxmpt];
double smat[MM][MM], mat[MM][MM];
double fped[NN], fibd[NN], fvr1[NN], fvr2[NN], fyang1[NN], fyang2[NN], fLH1[NN], fLH2[NN], fhom[NN];
double fibdBP[NN], fvr1BP[NN], fvr2BP[NN], fyang1BP[NN], fyang2BP[NN], fLH1BP[NN], fLH2BP[NN], fhomBP[NN];
double FITN[NN], FITN_F[NN], PHEN[NN], PHEN_F[NN];
double mean_FITN, mean_FITN_F, mean_PHEN, mean_PHEN_F;
double mean_FITN_0, mean_FITN_F_0, mean_PHEN_0, mean_PHEN_F_0;

double EHom, EHet, Hom[NN], sum_VR1[NN], sum_VR2[NN], sum_YANG1[NN], sum_YANG2[NN], sum_LH2[NN];
double EHomBP, EHetBP, HomBP[NN], sum_VR1BP[NN], sum_VR2BP[NN], sum_YANG1BP[NN], sum_YANG2BP[NN], sum_LH2BP[NN];

double sum_G, sum_F, sum_F2, sum_FG, var_F;

double NeFped, NeFibd, NeFvr1, NeFvr2, NeFyang1, NeFyang2, NeFLH1, NeFLH2, NeFhom;
struct acc AVE_NeFped, AVE_NeFibd, AVE_NeFvr1, AVE_NeFvr2, AVE_NeFyang1, AVE_NeFyang2, AVE_NeFLH1, AVE_NeFLH2, AVE_NeFhom;
double NeFpedBP, NeFibdBP, NeFvr1BP, NeFvr2BP, NeFyang1BP, NeFyang2BP, NeFLH1BP, NeFLH2BP, NeFhomBP;
struct acc AVE_NeFpedBP, AVE_NeFibdBP, AVE_NeFvr1BP, AVE_NeFvr2BP, AVE_NeFyang1BP, AVE_NeFyang2BP, AVE_NeFLH1BP, AVE_NeFLH2BP, AVE_NeFhomBP;

struct acc gmean_s[NN], gvar_s[NN], gmean_a[NN], gvar_a[NN], AVE_VA[NN], AVE_VD[NN], AVE_ID[NN], AVE_P[NN], AVE_VP[NN];
struct acc AVE_MAF[NN], AVE_D2[NN], AVE_r2[NN];
struct acc Fped[NN], Fibd[NN], Fvr1[NN], Fvr2[NN], Fyang1[NN], Fyang2[NN], FLH1[NN], FLH2[NN], Fhom[NN];
struct acc VFped[NN], VFibd[NN], VFvr1[NN], VFvr2[NN], VFyang1[NN], VFyang2[NN], VFLH1[NN], VFLH2[NN], VFhom[NN];
struct acc AVE_Fped[NN], AVE_Fibd[NN], AVE_Fvr1[NN], AVE_Fvr2[NN], AVE_Fyang1[NN], AVE_Fyang2[NN], AVE_FLH1[NN], AVE_FLH2[NN], AVE_Fhom[NN];
struct acc AVE_VFped[NN], AVE_VFibd[NN], AVE_VFvr1[NN], AVE_VFvr2[NN], AVE_VFyang1[NN], AVE_VFyang2[NN], AVE_VFLH1[NN], AVE_VFLH2[NN], AVE_VFhom[NN];

struct acc FpedBP[NN], FibdBP[NN], Fvr1BP[NN], Fvr2BP[NN], Fyang1BP[NN], Fyang2BP[NN], FLH1BP[NN], FLH2BP[NN], FhomBP[NN];
struct acc AVE_FpedBP[NN], AVE_FibdBP[NN], AVE_Fvr1BP[NN], AVE_Fvr2BP[NN], AVE_Fyang1BP[NN], AVE_Fyang2BP[NN], AVE_FLH1BP[NN], AVE_FLH2BP[NN], AVE_FhomBP[NN];
struct acc VFpedBP[NN], VFibdBP[NN], VFvr1BP[NN], VFvr2BP[NN], VFyang1BP[NN], VFyang2BP[NN], VFLH1BP[NN], VFLH2BP[NN], VFhomBP[NN];
struct acc AVE_VFpedBP[NN], AVE_VFibdBP[NN], AVE_VFvr1BP[NN], AVE_VFvr2BP[NN], AVE_VFyang1BP[NN], AVE_VFyang2BP[NN], AVE_VFLH1BP[NN], AVE_VFLH2BP[NN], AVE_VFhomBP[NN];

struct acc IDFped[NN], IDFibd[NN], IDFvr1[NN], IDFvr2[NN], IDFyang1[NN], IDFyang2[NN], IDFLH1[NN], IDFLH2[NN], IDFhom[NN];
struct acc AVE_IDFped, AVE_IDFibd, AVE_IDFvr1, AVE_IDFvr2, AVE_IDFyang1, AVE_IDFyang2, AVE_IDFLH1, AVE_IDFLH2, AVE_IDFhom;
struct acc IDFpedBP[NN], IDFibdBP[NN], IDFvr1BP[NN], IDFvr2BP[NN], IDFyang1BP[NN], IDFyang2BP[NN], IDFLH1BP[NN], IDFLH2BP[NN], IDFhomBP[NN];
struct acc AVE_IDFpedBP, AVE_IDFibdBP, AVE_IDFvr1BP, AVE_IDFvr2BP, AVE_IDFyang1BP, AVE_IDFyang2BP, AVE_IDFLH1BP, AVE_IDFLH2BP, AVE_IDFhomBP;

struct acc qi_0, qi_0_1, qi_1_2, qi_2_3, qi_3_4, qi_4_5, qi_5_6, qi_6_7, qi_7_8, qi_8_9, qi_9_10, qi_10;
struct acc AVE_qi_0, AVE_qi_0_1, AVE_qi_1_2, AVE_qi_2_3, AVE_qi_3_4, AVE_qi_4_5, AVE_qi_5_6, AVE_qi_6_7, AVE_qi_7_8, AVE_qi_8_9, AVE_qi_9_10, AVE_qi_10;

struct acc AVE_mean_FITN_0, AVE_mean_FITN_F_0, AVE_mean_PHEN_0, AVE_mean_PHEN_F_0;
struct acc AVE_mean_FITN, AVE_mean_FITN_F, AVE_mean_PHEN, AVE_mean_PHEN_F;
struct acc AVE_ID_FITN_0, AVE_ID_FITN, AVE_ID_PHEN_0, AVE_ID_PHEN;

struct acc AVE_EHom[NN], AVE_EHomBP[NN];
struct acc AVE_EHet[NN], AVE_EHetBP[NN];

struct covacc FpedFibd, FpedFvr1, FpedFvr2, FpedFyang1, FpedFyang2, FpedFLH1, FpedFLH2, FpedFhom;
struct covacc FibdFvr1, FibdFvr2, FibdFyang1, FibdFyang2, FibdFLH1, FibdFLH2, FibdFhom;
struct covacc Fvr1Fvr2, Fvr1Fyang1, Fvr1Fyang2, Fvr1FLH1, Fvr1FLH2, Fvr1Fhom;
struct covacc Fvr2Fyang1, Fvr2Fyang2, Fvr2FLH1, Fvr2FLH2, Fvr2Fhom;
struct covacc Fyang1Fyang2, Fyang1FLH1, Fyang1FLH2, Fyang1Fhom;
struct covacc Fyang2FLH1, Fyang2FLH2, Fyang2Fhom;
struct covacc FLH1FLH2, FLH1Fhom;
struct covacc FLH2Fhom;

struct acc AVE_FpedFibd, AVE_FpedFvr1, AVE_FpedFvr2, AVE_FpedFyang1, AVE_FpedFyang2, AVE_FpedFLH1, AVE_FpedFLH2, AVE_FpedFhom;
struct acc AVE_FibdFvr1, AVE_FibdFvr2, AVE_FibdFyang1, AVE_FibdFyang2, AVE_FibdFLH1, AVE_FibdFLH2, AVE_FibdFhom;
struct acc AVE_Fvr1Fvr2, AVE_Fvr1Fyang1, AVE_Fvr1Fyang2, AVE_Fvr1FLH1, AVE_Fvr1FLH2, AVE_Fvr1Fhom;
struct acc AVE_Fvr2Fyang1, AVE_Fvr2Fyang2, AVE_Fvr2FLH1, AVE_Fvr2FLH2, AVE_Fvr2Fhom;
struct acc AVE_Fyang1Fyang2, AVE_Fyang1FLH1, AVE_Fyang1FLH2, AVE_Fyang1Fhom;
struct acc AVE_Fyang2FLH1, AVE_Fyang2FLH2, AVE_Fyang2Fhom;
struct acc AVE_FLH1FLH2, AVE_FLH1Fhom;
struct acc AVE_FLH2Fhom;

struct covacc FpedFibdBP, FpedFvr1BP, FpedFvr2BP, FpedFyang1BP, FpedFyang2BP, FpedFLH1BP, FpedFLH2BP, FpedFhomBP;
struct covacc FibdFvr1BP, FibdFvr2BP, FibdFyang1BP, FibdFyang2BP, FibdFLH1BP, FibdFLH2BP, FibdFhomBP;
struct covacc Fvr1Fvr2BP, Fvr1Fyang1BP, Fvr1Fyang2BP, Fvr1FLH1BP, Fvr1FLH2BP, Fvr1FhomBP;
struct covacc Fvr2Fyang1BP, Fvr2Fyang2BP, Fvr2FLH1BP, Fvr2FLH2BP, Fvr2FhomBP;
struct covacc Fyang1Fyang2BP, Fyang1FLH1BP, Fyang1FLH2BP, Fyang1FhomBP;
struct covacc Fyang2FLH1BP, Fyang2FLH2BP, Fyang2FhomBP;
struct covacc FLH1FLH2BP, FLH1FhomBP;
struct covacc FLH2FhomBP;

struct acc AVE_FpedFibdBP, AVE_FpedFvr1BP, AVE_FpedFvr2BP, AVE_FpedFyang1BP, AVE_FpedFyang2BP, AVE_FpedFLH1BP, AVE_FpedFLH2BP, AVE_FpedFhomBP;
struct acc AVE_FibdFvr1BP, AVE_FibdFvr2BP, AVE_FibdFyang1BP, AVE_FibdFyang2BP, AVE_FibdFLH1BP, AVE_FibdFLH2BP, AVE_FibdFhomBP;
struct acc AVE_Fvr1Fvr2BP, AVE_Fvr1Fyang1BP, AVE_Fvr1Fyang2BP, AVE_Fvr1FLH1BP, AVE_Fvr1FLH2BP, AVE_Fvr1FhomBP;
struct acc AVE_Fvr2Fyang1BP, AVE_Fvr2Fyang2BP, AVE_Fvr2FLH1BP, AVE_Fvr2FLH2BP, AVE_Fvr2FhomBP;
struct acc AVE_Fyang1Fyang2BP, AVE_Fyang1FLH1BP, AVE_Fyang1FLH2BP, AVE_Fyang1FhomBP;
struct acc AVE_Fyang2FLH1BP, AVE_Fyang2FLH2BP, AVE_Fyang2FhomBP;
struct acc AVE_FLH1FLH2BP, AVE_FLH1FhomBP;
struct acc AVE_FLH2FhomBP;

FILE *fptr, *fOUTM, *fOUTV, *fgen, *frep, *fdat, *fpop, *ffmap, *ffreq, *ffped, *foutMF, *foutVF, *foutID, *foutNe, *foutR, *foutMFBP, *foutVFBP, *foutIDBP, *foutNeBP, *foutRBP, *foutINDF, *foutINDFBP, *fGRAPH;

/* ***************************************************** */

main()
{
	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");
	ffreq = fopen ("frequencies.dat","w");

	foutMF = fopen ("outfileF.dat","w");
	foutVF = fopen ("outfileVF.dat","w");
	foutID = fopen ("outfileID.dat","w");
	foutNe = fopen ("outfileNe.dat","w");
	foutR = fopen ("outfileR.dat","w");

	foutMFBP = fopen ("outfileFBP.dat","w");
	foutVFBP = fopen ("outfileVFBP.dat","w");
	foutIDBP = fopen ("outfileIDBP.dat","w");
	foutNeBP = fopen ("outfileNeBP.dat","w");
	foutRBP = fopen ("outfileRBP.dat","w");

	fOUTM = fopen ("outfileFLASTGEN.dat","w");
	fOUTV = fopen ("outfileVFLASTGEN.dat","w");
	foutINDF = fopen ("list_individual_Fs","w");
	foutINDFBP = fopen ("list_individual_FsBP","w");

	getinputs();
	headings();
	recombination_masks();

	/* POISSON RECOMBINATION NUMBER */
	if ( (exp(-(double)L) != 0.0)&&((double)L < normalthreshold) )
	generatepoissontable((double)L, &lastinrecpoissontable, recpoissontable, maxmpt-1);	

	natural_population();

	for (rep=1; rep<=replicates; rep++)
	{
		frep = fopen ("repfile.dat","a");
		fprintf (frep,"replicate %d\n", rep);
		fclose(frep);

		if (tracelevel!=0) fprintf (fptr,"\nreplicate %d\n\n", rep);

		sample();

		for (gen=0; gen<=generations; gen++)
		{
			frequency_genes ();
			linkage_disequilibrium ();
			genotypic_values ();

			phenotypeB ();
//			if (tracelevel!=0) dumpphenotypes();

			coancestry_matrix();
			if ((rep == 1) && (gen == generations))	plink_files();

			estimates_of_F();
			if (gen == 0)	ID_HOMOZYGOTES();

			if (gen == generations)
			{
				ID_HOMOZYGOTES();
				inbreeding_depression();
				correlations();
				estimates_Ne();
			}

			mating ();
			if (EC == 1)	disorder_parents ();
			else		selection ();

//			if (tracelevel!=0) dumpoffspring();
		}
		settozero();
	 	printout();
  	}
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();

	getintandskip("NINDNP (max 1000):",&NINDNP,2,1000);
	getintandskip("TNIND (max 1000):",&TNIND,2,1000);
	getintandskip("NIND (max 1000):",&NIND,2,1000);
	getrealandskip("sd_effect (normal with mean 0.0) :",&sd_effect,0.0,(double)infinity);
	getrealandskip("ha :",&ha,0.0,1.0);
	getrealandskip("VE :",&VE,0.0,(double)infinity);
	getrealandskip("Length of genome in Morgans (99:FreeRecom) :",&L,0.0,99.0);
	NLOCI=30;
	getintandskip("Number of generations :",&generations,0,200);
	getintandskip("Generation as BP :",&genBP,0,1000);
	getintandskip("Equalizing contributions (1) or not (0) :",&EC,0,1);
	getintandskip("Neutral for fitness (1) or not (0) :",&neutral,0,1);
	getintandskip("INITIALFREQ05 (yes 1 or not 0) :",&INITIALFREQ05,0,1);
	getrealandskip("MAF :",&MAF,0.0,1.0);
	getintandskip("Number of replicates :",&replicates,1,infinity);
}

/* **************************************************** */

headings()
{
	fgen = fopen ("genfile.dat","a");
	fprintf(fgen,"\nN=%d  N.LOCI=%d  reps=%d  gens=%d  L=%f\n", NIND, TOTLOCI, replicates, generations, L);
	fprintf(fgen,"\ngen   W+-SE             V(W)+-SE          q+-SE            HO+-SE           HS+-SE           HT+-SE           FST+-SE          NeH+-SE          NeFst+-SE\n");
	fclose(fgen);
}

/* ***************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ***************************************************** */

natural_population ()
{
	int x, numSNP, ds, g0, g1;
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
		posNP[k][l] = dpos;
		if (dps < -1.0) dps=(-1.0);
		if (da == -99.0) da=0.0;
		sNP[k][l] = dps;
		if (uniform() < 0.05)
		{
			atNP[k][l] = normal(0.0, sd_effect);
			if (atNP[k][l] < 0.0) atNP[k][l]=(-atNP[k][l]);
		}
		hsNP[k][l] = dh;
		hatNP[k][l] = ha;
//		if((k==4000)&&(l==0)) fprintf(fptr,"k=%d l=%d posNP=%llu chromNP=%d\n", k, l, posNP[k][l], chromNP[k][l]);
	}

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

/*	if (tracelevel!=0)	
	{
		fprintf(fptr, "\n");
		for (i=0; i<20; i++)
		{
			if ((gmNP[i][0][0] & RM[0])==RM[0])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
			if ((gmNP[i][0][1] & RM[0])==RM[0])	fprintf(fptr, "1 ");
			else								fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}
*/

	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));
		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);

		if (qNP[k][l] > 0.0)	NSEGLOCNP ++;
	} 

/*	if (tracelevel!=0)	
	{
		fprintf(fptr, "\n LEQ_NP = %f\n", LEQ_NP);
		for (k=0; k<NCRO; k++)	for (l=0; l<NLOCI; l++)
		fprintf(fptr, "\n k=%d l=%d   atNP=%f  hatNP=%f  qNP=%f", k, l, atNP[k][l], hatNP[k][l], qNP[k][l]);
	}
*/
	fclose(fdat);
}

/* ***************************************************** */

sample ()
{
	int g;

	/* ***** sample the first generation from the Base Population ***** */

	for (i=0; i<NINDNP; i++)	
	{
		ran_i = (int)(uniform() * NINDNP);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;

//			if (tracelevel!=0)	fprintf (fptr," i = %d  k = %d  gmNP0 = %d   gmNP1 = %d\n", i, k, gmNP[i][k][0], gmNP[i][k][1]);
		}
	}

	for (i=0; i<TNIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			gm[i][k][0]=gmNP[i][k][0];
			gm[i][k][1]=gmNP[i][k][1];

//			if (tracelevel!=0) for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0))	fprintf (fptr,"Sample %d    gm0 = %d   gm1 = %d\n", i, gm[i][k][0], gm[i][k][1]);
		}
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
//		if((k==4000)&&(l==0)) fprintf(fptr,"\n\nk=%d l=%d %llu chrom=%d\n\n", k, l, pos[k][l], chrom[k][l]);
		
//		if (tracelevel!=0)    if ((k<4)&&(l<5)) fprintf(fptr,"\n k=%d l=%d a=%f ha=%f", k, l, at[k][l], hat[k][l]);
	}
}

/* ***************************************************** */

frequency_genes ()
{
	double d_a, alfa_a;

	/* ***** MULTIALLELIC GENES ***** */

	if (gen == genBP)
	{
		for (k=0; k<NCRO; k++)
		for (i=0; i<TNIND; i++)
		{
			locus[i][k][0] = i;
			locus[i][k][1] = i+10000;

//			if (tracelevel!=0)	if (k<=3) fprintf (fptr,"locus=%d ind=%d   %d %d\n", k, i, locus[i][k][0], locus[i][k][1]);
		}
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<TNIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     	else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)TNIND)+(Aa/(2.0*(double)TNIND));
		if (gen == genBP)
		{
			qi[k][l] = (aa/(double)TNIND)+(Aa/(2.0*(double)TNIND));

			if (qi[k][l] == 0)							accum (&qi_0, 1.0);
			else if ((qi[k][l] > 0.0)&&(qi[k][l] <= 0.1))	accum (&qi_0_1, 1.0);
			else if ((qi[k][l] > 0.1)&&(qi[k][l] <= 0.2))	accum (&qi_1_2, 1.0);
			else if ((qi[k][l] > 0.2)&&(qi[k][l] <= 0.3))	accum (&qi_2_3, 1.0);
			else if ((qi[k][l] > 0.3)&&(qi[k][l] <= 0.4))	accum (&qi_3_4, 1.0);
			else if ((qi[k][l] > 0.4)&&(qi[k][l] <= 0.5))	accum (&qi_4_5, 1.0);
			else if ((qi[k][l] > 0.5)&&(qi[k][l] <= 0.6))	accum (&qi_5_6, 1.0);
			else if ((qi[k][l] > 0.6)&&(qi[k][l] <= 0.7))	accum (&qi_6_7, 1.0);
			else if ((qi[k][l] > 0.7)&&(qi[k][l] <= 0.8))	accum (&qi_7_8, 1.0);
			else if ((qi[k][l] > 0.8)&&(qi[k][l] <= 0.9))	accum (&qi_8_9, 1.0);
			else if ((qi[k][l] > 0.9)&&(qi[k][l] < 1.0))		accum (&qi_9_10, 1.0);
			else if (qi[k][l] == 1.0)					accum (&qi_10, 1.0);
		}

//		if (tracelevel!=0)    fprintf(fptr,"\n q[%d][%d]=%f", k, l, q[k][l]);

		if ((qi[k][l] != 0.0)&&(qi[k][l] != 1.0)&&(at[k][l] == 0.0)&&(s[k][l] == 0.0))
		{
			if (q[k][l] < (1.0-q[k][l]))	accum (&AVE_MAF[gen], q[k][l]);
			else						accum (&AVE_MAF[gen], 1.0-q[k][l]);
		}    
	}    

	if (gen == genBP)
	{
		fprintf(ffreq, "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
		accsum(&qi_0), accsum(&qi_0_1), accsum(&qi_1_2), accsum(&qi_2_3), accsum(&qi_3_4), accsum(&qi_4_5), accsum(&qi_5_6), accsum(&qi_6_7), accsum(&qi_7_8), accsum(&qi_8_9), accsum(&qi_9_10), accsum(&qi_10));
	
		accum (&AVE_qi_0, accsum(&qi_0));
		accum (&AVE_qi_0_1, accsum(&qi_0_1));
		accum (&AVE_qi_1_2, accsum(&qi_1_2));
		accum (&AVE_qi_2_3, accsum(&qi_2_3));
		accum (&AVE_qi_3_4, accsum(&qi_3_4));
		accum (&AVE_qi_4_5, accsum(&qi_4_5));
		accum (&AVE_qi_5_6, accsum(&qi_5_6));
		accum (&AVE_qi_6_7, accsum(&qi_6_7));
		accum (&AVE_qi_7_8, accsum(&qi_7_8));
		accum (&AVE_qi_8_9, accsum(&qi_8_9));
		accum (&AVE_qi_9_10, accsum(&qi_9_10));
		accum (&AVE_qi_10, accsum(&qi_10));
	}

	/* ******************* ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE ***************** */

	if (neutral == 1)
	{
		va = 0.0;
		vd = 0.0;
		id = 0.0;
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(at[k][l] != 0.0))
		{
			d_a = (at[k][l]/2.0) * (2.0*hat[k][l] - 1.0);
			if (at[k][l] >= 0.0)	alfa_a = (at[k][l]/2.0) + ( d_a * (1.0 - 2.0*q[k][l]) );
			else			alfa_a = (-at[k][l]/2.0) + ( d_a * (2.0*q[k][l] - 1.0) );

			va += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a * alfa_a;
			vd += (2.0 * d_a * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a * q[k][l] * (1.0-q[k][l]));
			id += (2.0 * d_a * q[k][l] * (1.0-q[k][l]));
		}    
//		if (tracelevel!=0)    fprintf(fptr,"\n va=%f vd=%f id=%f\n", va, vd, id);
		accum (&AVE_VA[gen], va);
		accum (&AVE_VD[gen], vd);
		accum (&AVE_ID[gen], id);
		if (gen == genBP)	idgenBP = id;
		if (gen == generations)	idgenLAST = id;
	}
	else
	{
		va = 0.0;
		vd = 0.0;
		id = 0.0;
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
		{
			d_a = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			if (s[k][l] >= 0.0)	alfa_a = (s[k][l]/2.0) + ( d_a * (1.0 - 2.0*q[k][l]) );
			else			alfa_a = (-s[k][l]/2.0) + ( d_a * (2.0*q[k][l] - 1.0) );

			va += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a * alfa_a;
			vd += (2.0 * d_a * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a * q[k][l] * (1.0-q[k][l]));
			id += (2.0 * d_a * q[k][l] * (1.0-q[k][l]));
		}    
//		if (tracelevel!=0)    fprintf(fptr,"\n va=%f vd=%f id=%f\n", va, vd, id);
		accum (&AVE_VA[gen], va);
		accum (&AVE_VD[gen], vd);
		accum (&AVE_ID[gen], id);
		if (gen == genBP)	idgenBP = id;
		if (gen == generations)	idgenLAST = id;
	}
}

/* ***************************************************** */

linkage_disequilibrium ()
{
	double D, VV, r2;
	double gam11, gam12, gam21, gam22; 

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI-1; l+=2)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(q[k][l+1] != 0.0)&&(q[k][l+1] != 1.0)&&(at[k][l] == 0.0)&&(at[k][l+1] == 0.0)&&(s[k][l] == 0.0)&&(s[k][l+1] == 0.0))
	{
		gam11=0.0, gam12=0.0, gam21=0.0, gam22=0.0;

		for (i=0; i<TNIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][0] & RM[l+1])==RM[l+1]))		gam11 ++;
			else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][0] & RM[l+1])!=RM[l+1]))	gam12 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][0] & RM[l+1])==RM[l+1]))	gam21 ++;
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][0] & RM[l+1])!=RM[l+1]))	gam22 ++;

			if (((gm[i][k][1] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l+1])==RM[l+1]))		gam11 ++;
			else if (((gm[i][k][1] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l+1])!=RM[l+1]))	gam12 ++;
			else if (((gm[i][k][1] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l+1])==RM[l+1]))	gam21 ++;
			else if (((gm[i][k][1] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l+1])!=RM[l+1]))	gam22 ++;
		}    

		D = ( (gam11/(2.0*TNIND)) * (gam22/(2.0*TNIND)) ) - ( (gam12/(2.0*TNIND)) * (gam21/(2.0*TNIND)) );
		VV = (q[k][l]*(1.0-q[k][l])) * (q[k][l+1]*(1.0-q[k][l+1]));
		if (VV != 0.0)	r2 = D*D / VV;
		if (VV != 0.0)	accum (&AVE_D2[gen], D*D);
		accum (&AVE_r2[gen], r2);
	}
}

/* ***************************************************** */

genotypic_values ()
{
//	if (tracelevel!=0)	fprintf(fptr,"\n Genotypic values \n");

	for (i=0; i<TNIND; i++)
	{
		pm_s[i] = 1.0;
		geno_a[i] = 0.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
//			if (tracelevel!=0)  if (i==12)  fprintf(fptr,"\n k=%d l=%d s=%f a=%f hs=%f ha=%f initialgen=%d",
//				k, l, s[k][l], at[k][l], hs[k][l], hat[k][l], initialgen[k][l]);

	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= (1.0 + s[k][l]);
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));

	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	geno_a[i] += at[k][l];
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else	geno_a[i] += (at[k][l]*hat[k][l]);
		}
		if (tracelevel!=0)
		{
			fprintf(fptr,"genotypic_values %d    pm_s = %f    geno_a = %f    pm_a = %f\n", i, pm_s[i], geno_a[i], pm_a[i]);
			for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0)) fprintf (fptr,"%d   gm0=%d   gm1=%d\n", i, gm[i][k][0], gm[i][k][1]);
		}
	}
}

/* ***************************************************** */

phenotypeB ()
{
	int ii, it;
	double gsum_s=0.0, gsum2_s=0.0;
	double gsum_a=0.0, gsum2_a=0.0;
	double gsum_pm_a=0.0, gsum2_pm_a=0.0;

	if (gen==0)	for (i=0; i<TNIND; i++)	pm_a[i] = geno_a[i] + normal(0.0, sqrt(VE));

	for (i=0; i<TNIND; i++)
	{
		gsum_s += pm_s[i];
		gsum2_s += (pm_s[i]*pm_s[i]);
		gsum_a += geno_a[i];
		gsum2_a += (geno_a[i]*geno_a[i]);
		gsum_pm_a += pm_a[i];
		gsum2_pm_a += (pm_a[i]*pm_a[i]);
	}

	accum (&gmean_s[gen], gsum_s/(double)TNIND);
	accum (&gvar_s[gen], (gsum2_s - (gsum_s*gsum_s / (double)TNIND)) / ((double)TNIND - 1.0));
	accum (&gmean_a[gen], gsum_a/(double)TNIND);
	accum (&gvar_a[gen], (gsum2_a - (gsum_a*gsum_a / (double)TNIND)) / ((double)TNIND - 1.0));
	accum (&AVE_P[gen], gsum_pm_a/(double)TNIND);
	accum (&AVE_VP[gen], (gsum2_pm_a - (gsum_pm_a*gsum_pm_a / (double)TNIND)) / ((double)TNIND - 1.0));

//	if (tracelevel!=0)   fprintf(fptr,"\ngmean_s = %f  gvar_s = %f C2 = %f gmean_a = %f  gvar_a = %f \n",
//	gsum_s/(double)TNIND, (gsum2_s - (gsum_s*gsum_s / (double)TNIND)) / (double)TNIND, ( (gsum2_s*(double)TNIND) / (gsum_s*gsum_s) ) - 1.0,
//	gsum_a/(double)TNIND, (gsum2_a - (gsum_a*gsum_a / (double)TNIND)) / (double)TNIND);
}

/* ***************************************************** */

dumpphenotypes()
{
	if (tracelevel==0)   return (0);

//	fprintf(fptr,"\n Fitness, Genotypic and Phenotypic values\n");
//	for (i=0; i<TNIND; i++)   fprintf(fptr,"i=%d pm_s=%f geno_a=%f pm_a=%f\n", i, pm_s[i], geno_a[i], pm_a[i]);
}

/* ***************************************************** */

plink_files()
{
	int last, lastpos;

	ffped = fopen ("data.ped","w");
	ffmap = fopen ("data.map","w");

	// data.map

	last = 1;
	lastpos = 1;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(at[k][l]==0.0))
	{
		if (chrom[k][l] != last)
		{
			last = chrom[k][l]; 
			lastpos = pos[k][l]; 
		}
		fprintf(ffmap,"%d SNP%llu 0 %llu\n", chrom[k][l], pos[k][l], pos[k][l]-lastpos);
	}

	for (i=0; i<TNIND; i++)
	{
		// data.ped

		fprintf(ffped,"1 IND%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(at[k][l]==0.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffped,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffped,"A T ");
		}
		fprintf(ffped,"\n");
	}
}

/* ************************************************************ */

coancestry_matrix()
{
	if (gen == genBP)
	{
		for (i=0; i<TNIND; i++)
		for (j=0; j<TNIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5;
			else		mat[i][j] = 0.0;
		}
	}
	else
	{
		for (i=0; i<TNIND; i++)
		for (j=0; j<TNIND; j++)
		{
			if (i == j)	mat[i][j] = 0.5 * (1.0 + smat[father[i]][mother[i]]);
			else		mat[i][j] = 0.25 * (smat[father[i]][father[j]] + smat[father[i]][mother[j]] + smat[mother[i]][father[j]] + smat[mother[i]][mother[j]]);
		}
	}

	for (i=0; i<TNIND; i++)
	for (j=0; j<TNIND; j++)
	{
		smat[i][j] = mat[i][j];
		if (j==0)
		{
			accum(&Fped[gen], 2.0*mat[i][i]-1.0); 
			accum(&FpedBP[gen], 2.0*mat[i][i]-1.0); 
			fped[i] = 2.0*mat[i][i]-1.0; 
		}
	}

/*	if (tracelevel!=0)    
	for (i=0; i<TNIND; i++)
	{
		fprintf(fptr, "gen=%d   i=%d   fat=%d   mot=%d   Fped=%f\n", gen, i, father[i], mother[i], 2.0*mat[i][i]-1.0);
	}
*/
}

/* ***************************************************** */

estimates_of_F ()
{
	// FREQUENCIES FROM CURRENT GENERATION

	numLOCI[gen] = 0.0;
	EHom = 0.0;
	EHet = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(at[k][l]==0.0)&&(q[k][l] >= MAF))
	{
		numLOCI[gen] ++;
		if (INITIALFREQ05 == 1)
		{
			EHom += 0.5;
			EHet += 0.5;
		}
		else
		{
			EHom += (1.0 - 2.0*q[k][l]*(1.0-q[k][l]));
			EHet += 2.0*q[k][l]*(1.0-q[k][l]);
		}
	}

	accum(&AVE_EHom[gen], EHom/numLOCI[gen]);
	accum(&AVE_EHet[gen], EHet/numLOCI[gen]);

	// FREQUENCIES FROM GENERATION BP

	numLOCIBP[gen] = 0.0;
	EHomBP = 0.0;
	EHetBP = 0.0;

	if ((rep == 1) && (gen == generations))
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		if ((qi[k][l] != 0.0)&&(qi[k][l] != 1.0)&&(at[k][l]==0.0)&&(q[k][l] >= MAF))
		fprintf(ffreq, "%f %f\n", qi[k][l], q[k][l]);
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((qi[k][l] != 0.0)&&(qi[k][l] != 1.0)&&(at[k][l]==0.0)&&(q[k][l] >= MAF))
	{
		numLOCIBP[gen] ++;
		EHomBP += (1.0 - 2.0*qi[k][l]*(1.0-qi[k][l]));
		EHetBP += 2.0*qi[k][l]*(1.0-qi[k][l]);
	}

	accum(&AVE_EHomBP[gen], EHomBP/numLOCIBP[gen]);
	accum(&AVE_EHetBP[gen], EHetBP/numLOCIBP[gen]);

	for (i=0; i<TNIND; i++)
	{
		// MULTILLELIC GENES

		fibd[i] = 0.0;
		for (k=0; k<NCRO; k++)
		{
			if (locus[i][k][0] == locus[i][k][1])
			{
				fibd[i] ++;
//				if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    locus[i][k][0]=%d locus[i][k][1]=%d \n", k, locus[i][k][0], locus[i][k][1]);
			}
		}
		fibd[i] = fibd[i] / NCRO;


		// FREQUENCIES FROM CURRENT GENERATION

		Hom[i] = 0.0;
		sum_VR1[i] = 0.0;
		sum_VR2[i] = 0.0; 
		sum_YANG1[i] = 0.0; 
		sum_YANG2[i] = 0.0; 
		sum_LH2[i] = 0.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(at[k][l]==0.0)&&(q[k][l] >= MAF))
		{
//			if (tracelevel!=0) if((k==0)&&(l<5))   fprintf(fptr,"\n q[%d][%d]=%f   EHom=%f", k, l, q[k][l], EHom);

			if (INITIALFREQ05 == 1)
			{
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
				{
					Hom[i] ++;
					sum_VR1[i] += (2.0-(2.0*0.5))*(2.0-(2.0*0.5));
					sum_VR2[i] += ( ( (2.0-2.0*0.5)*(2.0-2.0*0.5) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
					sum_YANG1[i] += ( (4.0)-((1.0+2.0*0.5)*2.0)+(2.0*0.5*0.5) ); 
					sum_YANG2[i] += ( (4.0)-((1.0+2.0*0.5)*2.0)+(2.0*0.5*0.5) ) / (2.0*0.5*(1.0-0.5)); 
				}
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
				{
					Hom[i] ++;
					sum_VR1[i] += (0.0-(2.0*0.5))*(0.0-(2.0*0.5));
					sum_VR2[i] += ( ( (0.0-(2.0*0.5))*(0.0-(2.0*0.5)) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
					sum_YANG1[i] += (2.0*0.5*0.5); 
					sum_YANG2[i] += (2.0*0.5*0.5) / (2.0*0.5*(1.0-0.5)); 
				}
				else
				{
					sum_VR1[i] += (1.0-(2.0*0.5))*(1.0-(2.0*0.5));
					sum_VR2[i] += ( ( (1.0-(2.0*0.5))*(1.0-(2.0*0.5)) ) / (2.0*0.5*(1.0-0.5)) ) - 1.0; 
					sum_YANG1[i] += ( (1.0)-((1.0+2.0*0.5)*1.0)+(2.0*0.5*0.5) ); 
					sum_YANG2[i] += ( (1.0)-((1.0+2.0*0.5)*1.0)+(2.0*0.5*0.5) ) / (2.0*0.5*(1.0-0.5)); 
					sum_LH2[i] += ( 1.0 / (2.0*0.5*(1.0-0.5)) ); 
				}
			}
			else
			{
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
				{
					Hom[i] ++;
					sum_VR1[i] += (2.0-(2.0*q[k][l]))*(2.0-(2.0*q[k][l]));
					sum_VR2[i] += ( ( (2.0-2.0*q[k][l])*(2.0-2.0*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
					sum_YANG1[i] += ( (4.0)-((1.0+2.0*q[k][l])*2.0)+(2.0*q[k][l]*q[k][l]) ); 
					sum_YANG2[i] += ( (4.0)-((1.0+2.0*q[k][l])*2.0)+(2.0*q[k][l]*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l])); 
				}
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
				{
					Hom[i] ++;
					sum_VR1[i] += (0.0-(2.0*q[k][l]))*(0.0-(2.0*q[k][l]));
					sum_VR2[i] += ( ( (0.0-(2.0*q[k][l]))*(0.0-(2.0*q[k][l])) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
					sum_YANG1[i] += (2.0*q[k][l]*q[k][l]); 
					sum_YANG2[i] += (2.0*q[k][l]*q[k][l]) / (2.0*q[k][l]*(1.0-q[k][l])); 
				}
				else
				{
					sum_VR1[i] += (1.0-(2.0*q[k][l]))*(1.0-(2.0*q[k][l]));
					sum_VR2[i] += ( ( (1.0-(2.0*q[k][l]))*(1.0-(2.0*q[k][l])) ) / (2.0*q[k][l]*(1.0-q[k][l])) ) - 1.0; 
					sum_YANG1[i] += ( (1.0)-((1.0+2.0*q[k][l])*1.0)+(2.0*q[k][l]*q[k][l]) ); 
					sum_YANG2[i] += ( (1.0)-((1.0+2.0*q[k][l])*1.0)+(2.0*q[k][l]*q[k][l]) ) / (2.0*q[k][l]*(1.0-q[k][l])); 
					sum_LH2[i] += ( 1.0 / (2.0*q[k][l]*(1.0-q[k][l])) ); 
				}
			}
		}

		fvr1[i] = (sum_VR1[i]/EHet)-1.0;
		fvr2[i] = sum_VR2[i] / numLOCI[gen];
		fyang1[i] = sum_YANG1[i] / EHet;
		fyang2[i] = sum_YANG2[i] / numLOCI[gen];
		fLH1[i] = (Hom[i]-EHom) / (numLOCI[gen]-EHom);
		fLH2[i] = 1.0 - ( sum_LH2[i] / numLOCI[gen]);
		fhom[i] = Hom[i] / numLOCI[gen];

		accum(&Fibd[gen], fibd[i]);
		accum(&Fvr1[gen], fvr1[i]);
		accum(&Fvr2[gen], fvr2[i]);
		accum(&Fyang1[gen], fyang1[i]);
		accum(&Fyang2[gen], fyang2[i]);
		accum(&FLH1[gen], fLH1[i]);
		accum(&FLH2[gen], fLH2[i]);
		accum(&Fhom[gen], fhom[i]);

		if ((neutral == 1) && (gen == generations)) fprintf(foutINDF,"%d %d %f %f %f %f %f %f %f %f %f %f\n", gen, i, pm_a[i], 2.0*mat[i][i]-1.0, fibd[i], fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);
		if ((neutral == 0) && (gen == generations)) fprintf(foutINDF,"%d %d %f %f %f %f %f %f %f %f %f %f\n", gen, i, pm_s[i], 2.0*mat[i][i]-1.0, fibd[i], fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);

		// MULTILLELIC GENES

		fibdBP[i] = 0.0;
		for (k=0; k<NCRO; k++)
		{
			if (locus[i][k][0] == locus[i][k][1])
			{
				fibdBP[i] ++;
//				if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    locus[i][k][0]=%d locus[i][k][1]=%d \n", k, locus[i][k][0], locus[i][k][1]);
			}
//			if (tracelevel!=0)	if (k<=3) fprintf(fptr," i=%d    k=%d    locus[i][k][0]=%d locus[i][k][1]=%d \n", i, k, locus[i][k][0], locus[i][k][1]);
		}
		fibdBP[i] = fibdBP[i] / NCRO;

		// FREQUENCIES FROM GENERATION BP

		HomBP[i] = 0.0;
		sum_VR1BP[i] = 0.0;
		sum_VR2BP[i] = 0.0; 
		sum_YANG1BP[i] = 0.0; 
		sum_YANG2BP[i] = 0.0; 
		sum_LH2BP[i] = 0.0; 

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((qi[k][l] != 0.0)&&(qi[k][l] != 1.0)&&(at[k][l]==0.0)&&(q[k][l] >= MAF))
		{
//			if (tracelevel!=0) if((k==0)&&(l<5))   fprintf(fptr,"\n q[%d][%d]=%f   EHomBP=%f", k, l, qi[k][l], EHomBP);

			// All SNPs

	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				HomBP[i] ++;
				sum_VR1BP[i] += (2.0-(2.0*qi[k][l]))*(2.0-(2.0*qi[k][l]));
				sum_VR2BP[i] += ( ( (2.0-2.0*qi[k][l])*(2.0-2.0*qi[k][l]) ) / (2.0*qi[k][l]*(1.0-qi[k][l])) ) - 1.0; 
				sum_YANG1BP[i] += ( (4.0)-((1.0+2.0*qi[k][l])*2.0)+(2.0*qi[k][l]*qi[k][l]) ); 
				sum_YANG2BP[i] += ( (4.0)-((1.0+2.0*qi[k][l])*2.0)+(2.0*qi[k][l]*qi[k][l]) ) / (2.0*qi[k][l]*(1.0-qi[k][l])); 
			}
			else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				HomBP[i] ++;
				sum_VR1BP[i] += (0.0-(2.0*qi[k][l]))*(0.0-(2.0*qi[k][l]));
				sum_VR2BP[i] += ( ( (0.0-(2.0*qi[k][l]))*(0.0-(2.0*qi[k][l])) ) / (2.0*qi[k][l]*(1.0-qi[k][l])) ) - 1.0; 
				sum_YANG1BP[i] += (2.0*qi[k][l]*qi[k][l]); 
				sum_YANG2BP[i] += (2.0*qi[k][l]*qi[k][l]) / (2.0*qi[k][l]*(1.0-qi[k][l])); 
			}
			else
			{
				sum_VR1BP[i] += (1.0-(2.0*qi[k][l]))*(1.0-(2.0*qi[k][l]));
				sum_VR2BP[i] += ( ( (1.0-(2.0*qi[k][l]))*(1.0-(2.0*qi[k][l])) ) / (2.0*qi[k][l]*(1.0-qi[k][l])) ) - 1.0; 
				sum_YANG1BP[i] += ( (1.0)-((1.0+2.0*qi[k][l])*1.0)+(2.0*qi[k][l]*qi[k][l]) ); 
				sum_YANG2BP[i] += ( (1.0)-((1.0+2.0*qi[k][l])*1.0)+(2.0*qi[k][l]*qi[k][l]) ) / (2.0*qi[k][l]*(1.0-qi[k][l])); 
				sum_LH2BP[i] += ( 1.0 / (2.0*qi[k][l]*(1.0-qi[k][l])) ); 
			}
		}

		// All SNPs
		fvr1BP[i] = (sum_VR1BP[i]/EHetBP)-1.0;
		fvr2BP[i] = sum_VR2BP[i] / numLOCIBP[gen];
		fyang1BP[i] = sum_YANG1BP[i] / EHetBP;
		fyang2BP[i] = sum_YANG2BP[i] / numLOCIBP[gen];
		fLH1BP[i] = (HomBP[i]-EHomBP) / (numLOCIBP[gen]-EHomBP);
		fLH2BP[i] = 1.0 - ( sum_LH2BP[i] / numLOCIBP[gen]);
		fhomBP[i] = Hom[i] / numLOCIBP[gen];

		accum(&FibdBP[gen], fibdBP[i]);
		accum(&Fvr1BP[gen], fvr1BP[i]);
		accum(&Fvr2BP[gen], fvr2BP[i]);
		accum(&Fyang1BP[gen], fyang1BP[i]);
		accum(&Fyang2BP[gen], fyang2BP[i]);
		accum(&FLH1BP[gen], fLH1BP[i]);
		accum(&FLH2BP[gen], fLH2BP[i]);
		accum(&FhomBP[gen], fhomBP[i]);

/*		if (tracelevel!=0)	
		{
			fprintf(fptr,"\ni=%d  sumVR1=%f  sumVR2=%f  sumYang1=%f  sumYang2=%f  Hom=%f  Ehom=%f\n", i, sum_VR1[i], sum_VR2[i], sum_YANG1[i], sum_YANG2[i], Hom[i], EHom);
			fprintf(fptr,"i=%d  numLOCI=%d  fped=%f  fibd=%f  fvr1=%f  fvr2=%f  fyang1=%f  fyang2=%f  fLH1=%f  fLH2=%f  fhom=%f\n", i, numLOCI[gen], fped[i], fibd[i], fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);

			fprintf(fptr,"\ni=%d  sumVR1BP=%f  sumVR2BP=%f  sumYang1BP=%f  sumYang2BP=%f  Hom=%f  EhomBP=%f\n", i, sum_VR1BP[i], sum_VR2BP[i], sum_YANG1BP[i], sum_YANG2BP[i], Hom[i], EHomBP);
			fprintf(fptr,"i=%d  numLOCIBP=%d  fped=%f  fibdBP=%f  fvr1BP=%f  fvr2BP=%f  fyang1BP=%f  fyang2BP=%f  fLH1BP=%f  fLH2BP=%f  fhomBP=%f\n", i, numLOCI[gen], fped[i], fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
		}
*/
//		if (tracelevel!=0) if (gen == generations) fprintf(fptr,"i=%d  pm_s=%f  numLOCI=%d  fped=%f  fvr1=%f  fvr2=%f  fyang1=%f  fyang2=%f  fLH1=%f  fLH2=%f  fhom=%f\n", i, pm_s[i], numLOCI[gen], fped[i], fvr1[i], fvr2[i], fyang1[i], fyang2[i], fLH1[i], fLH2[i], fhom[i]);

//		if (tracelevel!=0) if (rep==34) fprintf(fptr,"gen=%d  i=%d  mother=%d  father=%d  pm_s=%f  fped=%f  fibd=%f  fvr1=%f\n", gen, i, mother[i], father[i], pm_s[i], fped[i], fibd[i], fvr1[i]);
//		if ((tracelevel!=0)&&(gen == generations)&&(fped[i]>0.35)) fprintf(fptr,"i=%d  mother=%d  father=%d  pm_s=%f  fped=%f  fibd=%f  fvr1=%f\n", i, mother[i], father[i], pm_s[i], fped[i], fibd[i], fvr1[i]);

		if ((neutral == 1) && (gen == generations)) fprintf(foutINDFBP,"%d %d %f %f %f %f %f %f %f %f %f %f\n", gen, i, pm_a[i], 2.0*mat[i][i]-1.0, fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
		if ((neutral == 0) && (gen == generations)) fprintf(foutINDFBP,"%d %d %f %f %f %f %f %f %f %f %f %f\n", gen, i, pm_s[i], 2.0*mat[i][i]-1.0, fibdBP[i], fvr1BP[i], fvr2BP[i], fyang1BP[i], fyang2BP[i], fLH1BP[i], fLH2BP[i], fhomBP[i]);
	}

	accum(&AVE_Fped[gen], accmean(&Fped[gen]));
	accum(&AVE_Fibd[gen], accmean(&Fibd[gen]));
	accum(&AVE_Fvr1[gen], accmean(&Fvr1[gen]));
	accum(&AVE_Fvr2[gen], accmean(&Fvr2[gen]));
	accum(&AVE_Fyang1[gen], accmean(&Fyang1[gen]));
	accum(&AVE_Fyang2[gen], accmean(&Fyang2[gen]));
	accum(&AVE_FLH1[gen], accmean(&FLH1[gen]));
	accum(&AVE_FLH2[gen], accmean(&FLH2[gen]));
	accum(&AVE_Fhom[gen], accmean(&Fhom[gen]));

	accum(&AVE_VFped[gen], variance(&Fped[gen]));
	accum(&AVE_VFibd[gen], variance(&Fibd[gen]));
	accum(&AVE_VFvr1[gen], variance(&Fvr1[gen]));
	accum(&AVE_VFvr2[gen], variance(&Fvr2[gen]));
	accum(&AVE_VFyang1[gen], variance(&Fyang1[gen]));
	accum(&AVE_VFyang2[gen], variance(&Fyang2[gen]));
	accum(&AVE_VFLH1[gen], variance(&FLH1[gen]));
	accum(&AVE_VFLH2[gen], variance(&FLH2[gen]));
	accum(&AVE_VFhom[gen], variance(&Fhom[gen]));

	accum(&AVE_FpedBP[gen], accmean(&FpedBP[gen]));
	accum(&AVE_FibdBP[gen], accmean(&FibdBP[gen]));
	accum(&AVE_Fvr1BP[gen], accmean(&Fvr1BP[gen]));
	accum(&AVE_Fvr2BP[gen], accmean(&Fvr2BP[gen]));
	accum(&AVE_Fyang1BP[gen], accmean(&Fyang1BP[gen]));
	accum(&AVE_Fyang2BP[gen], accmean(&Fyang2BP[gen]));
	accum(&AVE_FLH1BP[gen], accmean(&FLH1BP[gen]));
	accum(&AVE_FLH2BP[gen], accmean(&FLH2BP[gen]));
	accum(&AVE_FhomBP[gen], accmean(&FhomBP[gen]));

	accum(&AVE_VFpedBP[gen], variance(&FpedBP[gen]));
	accum(&AVE_VFibdBP[gen], variance(&FibdBP[gen]));
	accum(&AVE_VFvr1BP[gen], variance(&Fvr1BP[gen]));
	accum(&AVE_VFvr2BP[gen], variance(&Fvr2BP[gen]));
	accum(&AVE_VFyang1BP[gen], variance(&Fyang1BP[gen]));
	accum(&AVE_VFyang2BP[gen], variance(&Fyang2BP[gen]));
	accum(&AVE_VFLH1BP[gen], variance(&FLH1BP[gen]));
	accum(&AVE_VFLH2BP[gen], variance(&FLH2BP[gen]));
	accum(&AVE_VFhomBP[gen], variance(&FhomBP[gen]));

	if (gen == generations)
	{
		fprintf(foutMF, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&Fped[gen]), accmean(&Fibd[gen]), accmean(&Fvr1[gen]), accmean(&Fvr2[gen]), accmean(&Fyang1[gen]), accmean(&Fyang2[gen]), accmean(&FLH1[gen]), accmean(&FLH2[gen]), accmean(&Fhom[gen]));

		fprintf(foutVF, "%f  %f  %f  %f  %f  %f  %f  %f  %f\n",
		variance(&Fped[gen]), variance(&Fibd[gen]), variance(&Fvr1[gen]), variance(&Fvr2[gen]), variance(&Fyang1[gen]), variance(&Fyang2[gen]), variance(&FLH1[gen]), variance(&FLH2[gen]), variance(&Fhom[gen]));

		fprintf(foutMFBP, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&FpedBP[gen]), accmean(&FibdBP[gen]), accmean(&Fvr1BP[gen]), accmean(&Fvr2BP[gen]), accmean(&Fyang1BP[gen]), accmean(&Fyang2BP[gen]), accmean(&FLH1BP[gen]), accmean(&FLH2BP[gen]), accmean(&FhomBP[gen]));

		fprintf(foutVFBP, "%f  %f  %f  %f  %f  %f  %f  %f  %f\n",
		variance(&FpedBP[gen]), variance(&FibdBP[gen]), variance(&Fvr1BP[gen]), variance(&Fvr2BP[gen]), variance(&Fyang1BP[gen]), variance(&Fyang2BP[gen]), variance(&FLH1BP[gen]), variance(&FLH2BP[gen]), variance(&FhomBP[gen]));
	}
}

/* ***************************************************** */

ID_HOMOZYGOTES()
{
	int rnd;

	for (i=0; i<TNIND; i++)
	{
		FITN[i] = 1.0;
		FITN_F[i] = 1.0;
		PHEN[i] = 0.0;
		PHEN_F[i] = 0.0;

		if (uniform() < 0.5) rnd = 1;
		else				 rnd = 2;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);
				PHEN[i] += at[k][l]*hat[k][l];

				if (rnd == 1)
				{
					FITN_F[i] *= (1.0 + s[k][l]);
					PHEN_F[i] += at[k][l];
				}
				else		/*11*/;
			}
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);
				PHEN[i] += at[k][l]*hat[k][l];

				if (rnd == 1)	/*11*/;
				else
				{
					FITN_F[i] *= (1.0 + s[k][l]);
					PHEN_F[i] += at[k][l];
				}
			}
			else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]);
				PHEN[i] += at[k][l];

				FITN_F[i] *= (1.0 + s[k][l]);
				PHEN_F[i] += at[k][l];
			}
		}
	}

	if (gen == 0)
	{
		mean_FITN_0 = 0.0;
		mean_FITN_F_0 = 0.0;
		mean_PHEN_0 = 0.0;
		mean_PHEN_F_0 = 0.0;

		for (i=0; i<TNIND; i++)
		{
			mean_FITN_0 += FITN[i]/(double)TNIND;
			mean_FITN_F_0 += FITN_F[i]/(double)TNIND;
			mean_PHEN_0 += PHEN[i]/(double)TNIND;
			mean_PHEN_F_0 += PHEN_F[i]/(double)TNIND;
		}
		accum(&AVE_mean_FITN_0, mean_FITN_0);
		accum(&AVE_mean_FITN_F_0, mean_FITN_F_0);
		accum(&AVE_mean_PHEN_0, mean_PHEN_0);
		accum(&AVE_mean_PHEN_F_0, mean_PHEN_F_0);

//printf("\n %f \n", mean_PHEN_0);
//printf("%f \n", accmean(&AVE_mean_PHEN_0));

		accum(&AVE_ID_FITN_0, -log(mean_FITN_F_0/mean_FITN_0));
		accum(&AVE_ID_PHEN_0, mean_PHEN_F_0-mean_PHEN_0);
//printf("%f \n", accmean(&AVE_ID_PHEN_0));
	}
	if (gen == generations)
	{
		mean_FITN = 0.0;
		mean_FITN_F = 0.0;
		mean_PHEN = 0.0;
		mean_PHEN_F = 0.0;

		for (i=0; i<TNIND; i++)
		{
			mean_FITN += FITN[i]/(double)TNIND;
			mean_FITN_F += FITN_F[i]/(double)TNIND;
			mean_PHEN += PHEN[i]/(double)TNIND;
			mean_PHEN_F += PHEN_F[i]/(double)TNIND;
		}

		accum(&AVE_mean_FITN, mean_FITN);
		accum(&AVE_mean_FITN_F, mean_FITN_F);
		accum(&AVE_mean_PHEN, mean_PHEN);
		accum(&AVE_mean_PHEN_F, mean_PHEN_F);

		accum(&AVE_ID_FITN, -log(mean_FITN_F/mean_FITN));
		accum(&AVE_ID_PHEN, mean_PHEN_F-mean_PHEN);
	}

	return(0);
}

/* **************************************************************************** */

inbreeding_depression ()
{
	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fped[i];
		sum_F2 += (fped[i] * fped[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fped[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fped[i] * log(pm_s[i]));
		}

//		fprintf(fptr,"i=%d pm_s=%f pm_a=%f Fped=%f sum_G=%f sum_F=%f sum_F2=%f sum_FG=%f var_F=%f ID=%f\n",
//		i, pm_s[i], pm_a[i], fped[i], sum_G, sum_F, sum_F2, sum_FG, sum_F2 - (sum_F * sum_F / (double)TNIND), (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F); 
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);

//	fprintf(fptr,"var_F=%f ID=%f\n",
//	var_F, (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F); 

	if (var_F > 0.0) accum (&IDFped[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);
	if (var_F > 0.0) accum (&IDFpedBP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fibd ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fibd[i];
		sum_F2 += (fibd[i] * fibd[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fibd[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fibd[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFibd[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fvr1 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fvr1[i];
		sum_F2 += (fvr1[i] * fvr1[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fvr1[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fvr1[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFvr1[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fvr2 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fvr2[i];
		sum_F2 += (fvr2[i] * fvr2[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fvr2[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fvr2[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFvr2[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fyang1 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fyang1[i];
		sum_F2 += (fyang1[i] * fyang1[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fyang1[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fyang1[i] * log(pm_s[i]));
		}

		fprintf(fptr,"i=%d pm_s=%f pm_a=%f Fyang1=%f sum_G=%f sum_F=%f sum_F2=%f sum_FG=%f var_F=%f ID=%f\n",
		i, pm_s[i], pm_a[i], fyang1[i], sum_G, sum_F, sum_F2, sum_FG, sum_F2 - (sum_F * sum_F / (double)TNIND), (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F); 
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFyang1[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fyang2 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fyang2[i];
		sum_F2 += (fyang2[i] * fyang2[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fyang2[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fyang2[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFyang2[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FLH1 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fLH1[i];
		sum_F2 += (fLH1[i] * fLH1[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fLH1[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fLH1[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFLH1[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FLH2 ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fLH2[i];
		sum_F2 += (fLH2[i] * fLH2[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fLH2[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fLH2[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFLH2[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fhom ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fhom[i];
		sum_F2 += (fhom[i] * fhom[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fhom[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fhom[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFhom[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FibdBP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fibdBP[i];
		sum_F2 += (fibdBP[i] * fibdBP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fibdBP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fibdBP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFibdBP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fvr1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fvr1BP[i];
		sum_F2 += (fvr1BP[i] * fvr1BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fvr1BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fvr1BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFvr1BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fvr2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fvr2BP[i];
		sum_F2 += (fvr2BP[i] * fvr2BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fvr2BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fvr2BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFvr2BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fyang1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fyang1BP[i];
		sum_F2 += (fyang1BP[i] * fyang1BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fyang1BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fyang1BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFyang1BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** Fyang2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fyang2BP[i];
		sum_F2 += (fyang2BP[i] * fyang2BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fyang2BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fyang2BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFyang2BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FLH1BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fLH1BP[i];
		sum_F2 += (fLH1BP[i] * fLH1BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fLH1BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fLH1BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFLH1BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FLH2BP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fLH2BP[i];
		sum_F2 += (fLH2BP[i] * fLH2BP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fLH2BP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fLH2BP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFLH2BP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************** FhomBP ************************

	sum_G=0.0; sum_F=0.0; sum_F2=0.0; sum_FG=0.0;
	for (i=0; i<TNIND; i++)
	{
		sum_F += fhomBP[i];
		sum_F2 += (fhomBP[i] * fhomBP[i]);
		if (neutral == 1)
		{
			sum_G += pm_a[i];
			sum_FG += (fhomBP[i] * pm_a[i]);
		}
		else
		{
			sum_G += log(pm_s[i]);
			sum_FG += (fhomBP[i] * log(pm_s[i]));
		}
	}
	var_F = sum_F2 - (sum_F * sum_F / (double)TNIND);
	if (var_F > 0.0) accum (&IDFhomBP[gen], (sum_FG - (sum_G * sum_F / (double)TNIND)) / var_F);

	// *********************************************************

	accum(&AVE_IDFped, accmean(&IDFped[gen]));
	accum(&AVE_IDFibd, accmean(&IDFibd[gen]));
	accum(&AVE_IDFvr1, accmean(&IDFvr1[gen]));
	accum(&AVE_IDFvr2, accmean(&IDFvr2[gen]));
	accum(&AVE_IDFyang1, accmean(&IDFyang1[gen]));
	accum(&AVE_IDFyang2, accmean(&IDFyang2[gen]));
	accum(&AVE_IDFLH1, accmean(&IDFLH1[gen]));
	accum(&AVE_IDFLH2, accmean(&IDFLH2[gen]));
	accum(&AVE_IDFhom, accmean(&IDFhom[gen]));

	accum(&AVE_IDFpedBP, accmean(&IDFpedBP[gen]));
	accum(&AVE_IDFibdBP, accmean(&IDFibdBP[gen]));
	accum(&AVE_IDFvr1BP, accmean(&IDFvr1BP[gen]));
	accum(&AVE_IDFvr2BP, accmean(&IDFvr2BP[gen]));
	accum(&AVE_IDFyang1BP, accmean(&IDFyang1BP[gen]));
	accum(&AVE_IDFyang2BP, accmean(&IDFyang2BP[gen]));
	accum(&AVE_IDFLH1BP, accmean(&IDFLH1BP[gen]));
	accum(&AVE_IDFLH2BP, accmean(&IDFLH2BP[gen]));
	accum(&AVE_IDFhomBP, accmean(&IDFhomBP[gen]));

	if (neutral==1)
	{
		fprintf(foutID, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			-idgenLAST, log(mean_FITN_F/mean_FITN), (mean_PHEN_F-mean_PHEN), accmean(&IDFvr1[gen]), accmean(&IDFvr2[gen]), accmean(&IDFyang1[gen]), accmean(&IDFyang2[gen]), accmean(&IDFLH1[gen]), accmean(&IDFLH2[gen]), accmean(&IDFhom[gen]));
		fprintf(foutIDBP, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			-idgenBP, log(mean_FITN_F_0/mean_FITN_0), (mean_PHEN_F_0-mean_PHEN_0), accmean(&IDFpedBP[gen]), accmean(&IDFibdBP[gen]), accmean(&IDFvr1BP[gen]), accmean(&IDFvr2BP[gen]), accmean(&IDFyang1BP[gen]), accmean(&IDFyang2BP[gen]), accmean(&IDFLH1BP[gen]), accmean(&IDFLH2BP[gen]), accmean(&IDFhomBP[gen]));
	}
	else
	{
		fprintf(foutID, "%6.4f  %6.4f  0.0000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			-idgenLAST, log(mean_FITN_F/mean_FITN), accmean(&IDFvr1[gen]), accmean(&IDFvr2[gen]), accmean(&IDFyang1[gen]), accmean(&IDFyang2[gen]), accmean(&IDFLH1[gen]), accmean(&IDFLH2[gen]), accmean(&IDFhom[gen]));
		fprintf(foutIDBP, "%6.4f  %6.4f  0.0000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			-idgenBP, log(mean_FITN_F_0/mean_FITN_0), accmean(&IDFpedBP[gen]), accmean(&IDFibdBP[gen]), accmean(&IDFvr1BP[gen]), accmean(&IDFvr2BP[gen]), accmean(&IDFyang1BP[gen]), accmean(&IDFyang2BP[gen]), accmean(&IDFLH1BP[gen]), accmean(&IDFLH2BP[gen]), accmean(&IDFhomBP[gen]));
	}
}

/* ***************************************************** */

mating ()
{
	int a, p1, p2, EE[MM], FF[MM], last;
	int numberrecs, nr, pointrec[MM][31], ncrorec[MM], rndk, rndl;

	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		sm[i][k][0]=gm[i][k][0];
		sm[i][k][1]=gm[i][k][1];

		// MULTIALLELIC GENES
		slocus[i][k][0] = locus[i][k][0];
		slocus[i][k][1] = locus[i][k][1];
	}

//	if (tracelevel!=0)	fprintf(fptr,"\n Parents \n");

	for (i=0; i<TNIND; i++)
	{
		generahijo: /* ***** */;

		if (EC == 1)
		{
			/***** parents EC *****/
			if ((i%2) == 0)	p1 = i;
			else			p1 = i-1;
			p2 = p1 + 1;
		}
		else
		{
			/***** parents RC *****/
			p1 = (int)(uniform()*NIND);
			do { p2 = (int)(uniform()*NIND); }
			while (p2 == p1);
		}
		father[i] = p1;
		mother[i] = p2;

		/*************************/

//		if (tracelevel!=0)   fprintf (fptr,"%d\t%d\n", p1, p2);

		if(L==99.0)
		{	    /* ******************* Free recombination ******************* */

			// MULTIALLELIC GENES
			for (k=0; k<NCRO; k++)
			{
				if (uniform() < 0.5)	locus[i][k][0] = slocus[p1][k][0];
				else					locus[i][k][0] = slocus[p1][k][1];

				if (uniform() < 0.5)	locus[i][k][1] = slocus[p2][k][0];
				else					locus[i][k][1] = slocus[p2][k][1];
			}

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
				FF[k] = ~EE[k];
			   	gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
			}
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);

			for (k=0; k<NCRO; k++)
			{
			   	EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
			   	FF[k] = ~EE[k];
			   	gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
			}
//			if (tracelevel!=0)   fprintf (fptr,"i=%d EE[0]=%d EE[1]=%d EE[2]=%d sm00=%d sm01=%d sm10=%d sm11=%d sm20=%d sm21=%d g00=%d g10=%d g20=%d \n", i, 			EE[0], EE[1], EE[2], sm[p1][0][0], sm[p1][0][1], sm[p1][1][0], sm[p1][1][1], 				sm[p1][2][0], sm[p1][2][1], gm[i][0][0], gm[i][1][0], gm[i][2][0]);
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
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs); 

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
			}

/*			if (i == 0)
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (pointrec[k][l] == 1)
				if (tracelevel!=0)	fprintf (fptr,"pointrec k=%d l=%d\n", k, l);	
			}
*/
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
				gm[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));

				// MULTIALLELIC GENES
				if ((EE[k] & RM[0]) == RM[0])	locus[i][k][0] = slocus[p1][k][0];
				else						locus[i][k][0] = slocus[p1][k][1];

//				if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    i=%d    p1=%d    slocus[p1][k][0]=%d slocus[p1][k][1]=%d locus[i][k][0]=%d \n", k, i, p1, slocus[p1][k][0], slocus[p1][k][1], locus[i][k][0]);
			}
/*
			if (tracelevel!=0)
			if (i == 0)
			{
				fprintf (fptr,"EE\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((EE[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
					}
					fprintf (fptr,"\n");
				}
				fprintf (fptr,"FF\n");
				for (k=0; k<NCRO; k++)
				{
					for (l=0; l<NLOCI; l++)
					{
						if((FF[k]&RM[l])==RM[l])  	fprintf (fptr,"1 ");
						else			   	fprintf (fptr,"0 ");
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
//			if (tracelevel!=0)   fprintf (fptr,"numberrecs=%d\n",numberrecs); 

			for (nr=0; nr<numberrecs; nr++)
			{
				rndk = (int)(uniform()*NCRO);
				rndl = (int)(uniform()*NLOCI);
				ncrorec[rndk] = 1;
				pointrec[rndk][rndl] = 1;
//				if (tracelevel!=0)	fprintf (fptr,"Rec %d rndk=%d  rndl=%d\n", nr, rndk, rndl);
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
				gm[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));

				// MULTIALLELIC GENES
				if ((EE[k] & RM[0]) == RM[0])	locus[i][k][1] = slocus[p2][k][0];
				else						locus[i][k][1] = slocus[p2][k][1];

//				if (tracelevel!=0)	if (k<=3) fprintf(fptr," k=%d    i=%d    p2=%d    slocus[p2][k][0]=%d slocus[p2][k][1]=%d locus[i][k][1]=%d \n", k, i, p2, slocus[p2][k][0], slocus[p2][k][1], locus[i][k][1]);
//				if (tracelevel!=0)	if (k<=3) if ((EE[k] & RM[0]) == RM[0]) fprintf(fptr,"EE[%d]=1\n", k);
			}
		}

		/* *****Genotypic and Phenotypic value of offspring***** */

		geno_a[i] = 0.0;
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	geno_a[i] += at[k][l];
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 		/* AA */;
			else	geno_a[i] += (at[k][l]*hat[k][l]);
		}
		pm_a[i] = geno_a[i] + normal(0.0, sqrt(VE));
// 		if (tracelevel!=0)	fprintf (fptr," %d   geno_a = %f   pm_a = %f\n", i, geno_a[i], pm_a[i]);

		/* *****Fitness value of offspring***** */

    		if (neutral == 0)
		{
			pm_s[i] = 1.0;

			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= (1.0 + s[k][l]);
				else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 		/* AA */;
				else	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l]));
			}
		}
		else	pm_s[i] = 1.0;
// 		if (tracelevel!=0)	fprintf (fptr,"Offspring i = %d   pm_s = %f\n", i, pm_s[i]);

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

disorder_parents ()
{
	int a;
	double w;

	for (i=0; i<TNIND; i++)
	{
		ran_i=(int)(uniform()*TNIND);

		w=geno_a[i]; geno_a[i]=geno_a[ran_i]; geno_a[ran_i]=w;
		w=pm_a[i]; pm_a[i]=pm_a[ran_i]; pm_a[ran_i]=w;
		w=pm_s[i]; pm_s[i]=pm_s[ran_i]; pm_s[ran_i]=w;
		a=mother[i]; mother[i]=mother[ran_i]; mother[ran_i]=a;
		a=father[i]; father[i]=father[ran_i]; father[ran_i]=a;

		for (k=0; k<NCRO; k++)
		{
			a=gm[i][k][0]; gm[i][k][0]=gm[ran_i][k][0]; gm[ran_i][k][0]=a;
			a=gm[i][k][1]; gm[i][k][1]=gm[ran_i][k][1]; gm[ran_i][k][1]=a;

			// MULTIALLELIC GENES
			a=locus[i][k][0]; locus[i][k][0]=locus[ran_i][k][0]; locus[ran_i][k][0]=a;
			a=locus[i][k][1]; locus[i][k][1]=locus[ran_i][k][1]; locus[ran_i][k][1]=a;
		}
	}

/* 	if (tracelevel!=0)
	for (i=0; i<TNIND; i++)
	{
		fprintf (fptr,"Disordered parents %d   geno_a = %f   pm_a = %f\n", i, geno_a[i], pm_a[i]);
		for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0)) fprintf (fptr,"%d   gm0=%d   gm1=%d\n", i, gm[i][k][0], gm[i][k][1]);
	}
*/
}

/* ***************************************************** */

selection ()
{
	int a;
	double w;

	/* ******************* Order by phenotype ******************* */

	for (i=0; i<TNIND-1; i++)
	for (j=i+1; j<TNIND; j++)
	{
		if (pm_a[j] > pm_a[i])
		{
			w=geno_a[j]; geno_a[j]=geno_a[i]; geno_a[i]=w;
			w=pm_a[j]; pm_a[j]=pm_a[i]; pm_a[i]=w;
			w=pm_s[j]; pm_s[j]=pm_s[i]; pm_s[i]=w;
			a=mother[j]; mother[j]=mother[i]; mother[i]=a;
			a=father[j]; father[j]=father[i]; father[i]=a;

			for (k=0; k<NCRO; k++)
			{
				a=gm[j][k][0]; gm[j][k][0]=gm[i][k][0]; gm[i][k][0]=a;
				a=gm[j][k][1]; gm[j][k][1]=gm[i][k][1]; gm[i][k][1]=a;

				// MULTIALLELIC GENES
				a=locus[j][k][0]; locus[j][k][0]=locus[i][k][0]; locus[i][k][0]=a;
				a=locus[j][k][1]; locus[j][k][1]=locus[i][k][1]; locus[i][k][1]=a;
			}
		}
	}

/* 	if (tracelevel!=0)
	for (i=0; i<TNIND; i++)
	{
		fprintf (fptr,"Selected %d   geno_a = %f   pm_a = %f\n", i, geno_a[i], pm_a[i]);
		for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0)) fprintf (fptr,"%d   gm0=%d   gm1=%d\n", i, gm[i][k][0], gm[i][k][1]);
	}
*/
}

/* ***************************************************** */

correlations ()
{
	for (i=0; i<TNIND; i++)
	{
		covaccum (&FpedFibd, fped[i], fibd[i]);
		covaccum (&FpedFvr1, fped[i], fvr1[i]);
		covaccum (&FpedFvr2, fped[i], fvr2[i]);
		covaccum (&FpedFyang1, fped[i], fyang1[i]);
		covaccum (&FpedFyang2, fped[i], fyang2[i]);
		covaccum (&FpedFLH1, fped[i], fLH1[i]);
		covaccum (&FpedFLH2, fped[i], fLH2[i]);
		covaccum (&FpedFhom, fped[i], fhom[i]);
		covaccum (&FibdFvr1, fibd[i], fvr1[i]);
		covaccum (&FibdFvr2, fibd[i], fvr2[i]);
		covaccum (&FibdFyang1, fibd[i], fyang1[i]);
		covaccum (&FibdFyang2, fibd[i], fyang2[i]);
		covaccum (&FibdFLH1, fibd[i], fLH1[i]);
		covaccum (&FibdFLH2, fibd[i], fLH2[i]);
		covaccum (&FibdFhom, fibd[i], fhom[i]);
		covaccum (&Fvr1Fvr2, fvr1[i], fvr2[i]);
		covaccum (&Fvr1Fyang1, fvr1[i], fyang1[i]);
		covaccum (&Fvr1Fyang2, fvr1[i], fyang2[i]);
		covaccum (&Fvr1FLH1, fvr1[i], fLH1[i]);
		covaccum (&Fvr1FLH2, fvr1[i], fLH2[i]);
		covaccum (&Fvr1Fhom, fvr1[i], fhom[i]);
		covaccum (&Fvr2Fyang1, fvr2[i], fyang1[i]);
		covaccum (&Fvr2Fyang2, fvr2[i], fyang2[i]);
		covaccum (&Fvr2FLH1, fvr2[i], fLH1[i]);
		covaccum (&Fvr2FLH2, fvr2[i], fLH2[i]);
		covaccum (&Fvr2Fhom, fvr2[i], fhom[i]);
		covaccum (&Fyang1Fyang2, fyang1[i], fyang2[i]);
		covaccum (&Fyang1FLH1, fyang1[i], fLH1[i]);
		covaccum (&Fyang1FLH2, fyang1[i], fLH2[i]);
		covaccum (&Fyang1Fhom, fyang1[i], fhom[i]);
		covaccum (&Fyang2FLH1, fyang2[i], fLH1[i]);
		covaccum (&Fyang2FLH2, fyang2[i], fLH2[i]);
		covaccum (&Fyang2Fhom, fyang2[i], fhom[i]);
		covaccum (&FLH1FLH2, fLH1[i], fLH2[i]);
		covaccum (&FLH1Fhom, fLH1[i], fhom[i]);
		covaccum (&FLH2Fhom, fLH2[i], fhom[i]);

		covaccum (&FpedFibdBP, fped[i], fibdBP[i]);
		covaccum (&FpedFvr1BP, fped[i], fvr1BP[i]);
		covaccum (&FpedFvr2BP, fped[i], fvr2BP[i]);
		covaccum (&FpedFyang1BP, fped[i], fyang1BP[i]);
		covaccum (&FpedFyang2BP, fped[i], fyang2BP[i]);
		covaccum (&FpedFLH1BP, fped[i], fLH1BP[i]);
		covaccum (&FpedFLH2BP, fped[i], fLH2BP[i]);
		covaccum (&FpedFhomBP, fped[i], fhomBP[i]);
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

	accum (&AVE_FpedFibd, correlation(&FpedFibd));
	accum (&AVE_FpedFvr1, correlation(&FpedFvr1));
	accum (&AVE_FpedFvr2, correlation(&FpedFvr2));
	accum (&AVE_FpedFyang1, correlation(&FpedFyang1));
	accum (&AVE_FpedFyang2, correlation(&FpedFyang2));
	accum (&AVE_FpedFLH1, correlation(&FpedFLH1));
	accum (&AVE_FpedFLH2, correlation(&FpedFLH2));
	accum (&AVE_FpedFhom, correlation(&FpedFhom));
	accum (&AVE_Fvr1Fvr2, correlation(&Fvr1Fvr2));
	accum (&AVE_FibdFvr1, correlation(&FibdFvr1));
	accum (&AVE_FibdFvr2, correlation(&FibdFvr2));
	accum (&AVE_FibdFyang1, correlation(&FibdFyang1));
	accum (&AVE_FibdFyang2, correlation(&FibdFyang2));
	accum (&AVE_FibdFLH1, correlation(&FibdFLH1));
	accum (&AVE_FibdFLH2, correlation(&FibdFLH2));
	accum (&AVE_FibdFhom, correlation(&FibdFhom));
	accum (&AVE_Fvr1Fyang1, correlation(&Fvr1Fyang1));
	accum (&AVE_Fvr1Fyang2, correlation(&Fvr1Fyang2));
	accum (&AVE_Fvr1FLH1, correlation(&Fvr1FLH1));
	accum (&AVE_Fvr1FLH2, correlation(&Fvr1FLH2));
	accum (&AVE_Fvr1Fhom, correlation(&Fvr1Fhom));
	accum (&AVE_Fvr2Fyang1, correlation(&Fvr2Fyang1));
	accum (&AVE_Fvr2Fyang2, correlation(&Fvr2Fyang2));
	accum (&AVE_Fvr2FLH1, correlation(&Fvr2FLH1));
	accum (&AVE_Fvr2FLH2, correlation(&Fvr2FLH2));
	accum (&AVE_Fvr2Fhom, correlation(&Fvr2Fhom));
	accum (&AVE_Fyang1Fyang2, correlation(&Fyang1Fyang2));
	accum (&AVE_Fyang1FLH1, correlation(&Fyang1FLH1));
	accum (&AVE_Fyang1FLH2, correlation(&Fyang1FLH2));
	accum (&AVE_Fyang1Fhom, correlation(&Fyang1Fhom));
	accum (&AVE_Fyang2FLH1, correlation(&Fyang2FLH1));
	accum (&AVE_Fyang2FLH2, correlation(&Fyang2FLH2));
	accum (&AVE_Fyang2Fhom, correlation(&Fyang2Fhom));
	accum (&AVE_FLH1FLH2, correlation(&FLH1FLH2));
	accum (&AVE_FLH1Fhom, correlation(&FLH1Fhom));
	accum (&AVE_FLH2Fhom, correlation(&FLH2Fhom));

	fprintf(foutR, "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
	correlation(&FpedFibd), correlation(&FpedFvr1), correlation(&FpedFvr2), correlation(&FpedFyang1), correlation(&FpedFyang2), correlation(&FpedFLH1), correlation(&FpedFLH2), correlation(&FpedFhom),
	correlation(&FibdFvr1), correlation(&FibdFvr2), correlation(&FibdFyang1), correlation(&FibdFyang2), correlation(&FibdFLH1), correlation(&FibdFLH2), correlation(&FibdFhom),
	correlation(&Fvr1Fvr2), correlation(&Fvr1Fyang1), correlation(&Fvr1Fyang2), correlation(&Fvr1FLH1), correlation(&Fvr1FLH2), correlation(&Fvr1Fhom), correlation(&Fvr2Fyang1), correlation(&Fvr2Fyang2),
	correlation(&Fvr2FLH1), correlation(&Fvr2FLH2), correlation(&Fvr2Fhom), correlation(&Fyang1Fyang2), correlation(&Fyang1FLH1), correlation(&Fyang1FLH2), correlation(&Fyang1Fhom), correlation(&Fyang2FLH1), correlation(&Fyang2FLH2), correlation(&Fyang2Fhom), correlation(&FLH1FLH2), correlation(&FLH1Fhom), correlation(&FLH2Fhom));

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

	fprintf(foutRBP, "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
	correlation(&FpedFibdBP), correlation(&FpedFvr1BP), correlation(&FpedFvr2BP), correlation(&FpedFyang1BP), correlation(&FpedFyang2BP), correlation(&FpedFLH1BP), correlation(&FpedFLH2BP), correlation(&FpedFhomBP),
	correlation(&FibdFvr1BP), correlation(&FibdFvr2BP), correlation(&FibdFyang1BP), correlation(&FibdFyang2BP), correlation(&FibdFLH1BP), correlation(&FibdFLH2BP), correlation(&FibdFhomBP),
	correlation(&Fvr1Fvr2BP), correlation(&Fvr1Fyang1BP), correlation(&Fvr1Fyang2BP), correlation(&Fvr1FLH1BP), correlation(&Fvr1FLH2BP), correlation(&Fvr1FhomBP), correlation(&Fvr2Fyang1BP), correlation(&Fvr2Fyang2BP),
	correlation(&Fvr2FLH1BP), correlation(&Fvr2FLH2BP), correlation(&Fvr2FhomBP), correlation(&Fyang1Fyang2BP), correlation(&Fyang1FLH1BP), correlation(&Fyang1FLH2BP), correlation(&Fyang1FhomBP), correlation(&Fyang2FLH1BP), correlation(&Fyang2FLH2BP), correlation(&Fyang2FhomBP), correlation(&FLH1FLH2BP), correlation(&FLH1FhomBP), correlation(&FLH2FhomBP));
}

/* ***************************************************** */

estimates_Ne ()
{
	NeFpedBP = (5.0) / ( 2.0*(accmean(&FpedBP[generations])-accmean(&FpedBP[generations-5]))/(1.0-accmean(&FpedBP[generations-5])) );
	NeFibdBP = (5.0) / ( 2.0*(accmean(&FibdBP[generations])-accmean(&FibdBP[generations-5]))/(1.0-accmean(&FibdBP[generations-5])) );
	NeFvr1BP = (5.0) / ( 2.0*(accmean(&Fvr1BP[generations])-accmean(&Fvr1BP[generations-5]))/(1.0-accmean(&Fvr1BP[generations-5])) );
	NeFvr2BP = (5.0) / ( 2.0*(accmean(&Fvr2BP[generations])-accmean(&Fvr2BP[generations-5]))/(1.0-accmean(&Fvr2BP[generations-5])) );
	NeFyang1BP = (5.0) / ( 2.0*(accmean(&Fyang1BP[generations])-accmean(&Fyang1BP[generations-5]))/(1.0-accmean(&Fyang1BP[generations-5])) );
	NeFyang2BP = (5.0) / ( 2.0*(accmean(&Fyang2BP[generations])-accmean(&Fyang2BP[generations-5]))/(1.0-accmean(&Fyang2BP[generations-5])) );
	NeFLH1BP = (5.0) / ( 2.0*(accmean(&FLH1BP[generations])-accmean(&FLH1BP[generations-5]))/(1.0-accmean(&FLH1BP[generations-5])) );
	NeFLH2BP = (5.0) / ( 2.0*(accmean(&FLH2BP[generations])-accmean(&FLH2BP[generations-5]))/(1.0-accmean(&FLH2BP[generations-5])) );
	NeFhomBP = (5.0) / ( 2.0*(accmean(&FhomBP[generations])-accmean(&FhomBP[generations-5]))/(1.0-accmean(&FhomBP[generations-5])) );

	accum (&AVE_NeFpedBP, NeFpedBP);
	accum (&AVE_NeFibdBP, NeFibdBP);
	accum (&AVE_NeFvr1BP, NeFvr1BP);
	accum (&AVE_NeFvr2BP, NeFvr2BP);
	accum (&AVE_NeFyang1BP, NeFyang1BP);
	accum (&AVE_NeFyang2BP, NeFyang2BP);
	accum (&AVE_NeFLH1BP, NeFLH1BP);
	accum (&AVE_NeFLH2BP, NeFLH2BP);
	accum (&AVE_NeFhomBP, NeFhomBP);

	NeFped = (5.0) / ( 2.0*(accmean(&Fped[generations])-accmean(&Fped[generations-5]))/(1.0-accmean(&Fped[generations-5])) );
	NeFibd = (5.0) / ( 2.0*(accmean(&Fibd[generations])-accmean(&Fibd[generations-5]))/(1.0-accmean(&Fibd[generations-5])) );
	NeFvr1 = (5.0) / ( 2.0*(accmean(&Fvr1[generations])-accmean(&Fvr1[generations-5]))/(1.0-accmean(&Fvr1[generations-5])) );
	NeFvr2 = (5.0) / ( 2.0*(accmean(&Fvr2[generations])-accmean(&Fvr2[generations-5]))/(1.0-accmean(&Fvr2[generations-5])) );
	NeFyang1 = (5.0) / ( 2.0*(accmean(&Fyang1[generations])-accmean(&Fyang1[generations-5]))/(1.0-accmean(&Fyang1[generations-5])) );
	NeFyang2 = (5.0) / ( 2.0*(accmean(&Fyang2[generations])-accmean(&Fyang2[generations-5]))/(1.0-accmean(&Fyang2[generations-5])) );
	NeFLH1 = (5.0) / ( 2.0*(accmean(&FLH1[generations])-accmean(&FLH1[generations-5]))/(1.0-accmean(&FLH1[generations-5])) );
	NeFLH2 = (5.0) / ( 2.0*(accmean(&FLH2[generations])-accmean(&FLH2[generations-5]))/(1.0-accmean(&FLH2[generations-5])) );
	NeFhom = (5.0) / ( 2.0*(accmean(&Fhom[generations])-accmean(&Fhom[generations-5]))/(1.0-accmean(&Fhom[generations-5])) );

	accum (&AVE_NeFped, NeFped);
	accum (&AVE_NeFibd, NeFibd);
	accum (&AVE_NeFvr1, NeFvr1);
	accum (&AVE_NeFvr2, NeFvr2);
	accum (&AVE_NeFyang1, NeFyang1);
	accum (&AVE_NeFyang2, NeFyang2);
	accum (&AVE_NeFLH1, NeFLH1);
	accum (&AVE_NeFLH2, NeFLH2);
	accum (&AVE_NeFhom, NeFhom);

	fprintf(foutNe, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", NeFped, NeFibd, NeFvr1, NeFvr2, NeFyang1, NeFyang2, NeFLH1, NeFLH2, NeFhom);
	fprintf(foutNeBP, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", NeFpedBP, NeFibdBP, NeFvr1BP, NeFvr2BP, NeFyang1BP, NeFyang2BP, NeFLH1BP, NeFLH2BP, NeFhomBP);
}

/* ***************************************************** */

dumpoffspring()
{
	if (tracelevel==0)   return (0);

//	fprintf(fptr,"\n Offspring after selection (pm_a)\n");	
//	for (i=0; i<TNIND; i++)   fprintf(fptr,"%d pm_a=%f\n",i, pm_a[i]);
}

/* ***************************************************** */

settozero()
{
	for (gen=0; gen<=generations; gen++)
	{
		initacc (&Fped[gen]);
		initacc (&Fibd[gen]);
		initacc (&Fvr1[gen]);
		initacc (&Fvr2[gen]);
		initacc (&Fyang1[gen]);
		initacc (&Fyang2[gen]);
		initacc (&FLH1[gen]);
		initacc (&FLH2[gen]);
		initacc (&Fhom[gen]);

		initacc (&FpedBP[gen]);
		initacc (&FibdBP[gen]);
		initacc (&Fvr1BP[gen]);
		initacc (&Fvr2BP[gen]);
		initacc (&Fyang1BP[gen]);
		initacc (&Fyang2BP[gen]);
		initacc (&FLH1BP[gen]);
		initacc (&FLH2BP[gen]);
		initacc (&FhomBP[gen]);

		initacc (&IDFped[gen]);
		initacc (&IDFibd[gen]);
		initacc (&IDFvr1[gen]);
		initacc (&IDFvr2[gen]);
		initacc (&IDFyang1[gen]);
		initacc (&IDFyang2[gen]);
		initacc (&IDFLH1[gen]);
		initacc (&IDFLH2[gen]);
		initacc (&IDFhom[gen]);

		initacc (&IDFpedBP[gen]);
		initacc (&IDFibdBP[gen]);
		initacc (&IDFvr1BP[gen]);
		initacc (&IDFvr2BP[gen]);
		initacc (&IDFyang1BP[gen]);
		initacc (&IDFyang2BP[gen]);
		initacc (&IDFLH1BP[gen]);
		initacc (&IDFLH2BP[gen]);
		initacc (&IDFhomBP[gen]);
	}

	initacc (&qi_0);
	initacc (&qi_0_1);
	initacc (&qi_1_2);
	initacc (&qi_2_3);
	initacc (&qi_3_4);
	initacc (&qi_4_5);
	initacc (&qi_5_6);
	initacc (&qi_6_7);
	initacc (&qi_7_8);
	initacc (&qi_8_9);
	initacc (&qi_9_10);
	initacc (&qi_10);

	initcovacc (&FpedFibd);
	initcovacc (&FpedFvr1);
	initcovacc (&FpedFvr2);
	initcovacc (&FpedFyang1);
	initcovacc (&FpedFyang2);
	initcovacc (&FpedFLH1);
	initcovacc (&FpedFLH2);
	initcovacc (&FpedFhom);
	initcovacc (&FibdFvr1);
	initcovacc (&FibdFvr2);
	initcovacc (&FibdFyang1);
	initcovacc (&FibdFyang2);
	initcovacc (&FibdFLH1);
	initcovacc (&FibdFLH2);
	initcovacc (&FibdFhom);
	initcovacc (&Fvr1Fvr2);
	initcovacc (&Fvr1Fyang1);
	initcovacc (&Fvr1Fyang2);
	initcovacc (&Fvr1FLH1);
	initcovacc (&Fvr1FLH2);
	initcovacc (&Fvr1Fhom);
	initcovacc (&Fvr2Fyang1);
	initcovacc (&Fvr2Fyang2);
	initcovacc (&Fvr2FLH1);
	initcovacc (&Fvr2FLH2);
	initcovacc (&Fvr2Fhom);
	initcovacc (&Fyang1Fyang2);
	initcovacc (&Fyang1FLH1);
	initcovacc (&Fyang1FLH2);
	initcovacc (&Fyang1Fhom);
	initcovacc (&Fyang2FLH1);
	initcovacc (&Fyang2FLH2);
	initcovacc (&Fyang2Fhom);
	initcovacc (&FLH1FLH2);
	initcovacc (&FLH1Fhom);
	initcovacc (&FLH2Fhom);

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
	fgen = fopen ("genfile.dat","w");
	fGRAPH = fopen ("GRAPHS","w");

	fprintf(fgen, "TNIND=%d  Nsel=%d  L=%f  NSEGLOCNP=%d\n", TNIND, NIND, L, NSEGLOCNP);

	fprintf(fgen, "\n*********** Mean fitness and phenotype ***********\n\n");

	fprintf(fgen, "\ngen  W         varW      Phe       VP       VA         VD         ID\n");
	for(gen=0; gen<=generations; gen++)
		fprintf(fgen, "%3d  %f  %f  %f  %f  %f  %f  %f\n",
		gen, accmean(&gmean_s[gen]), accmean(&gvar_s[gen]), accmean(&AVE_P[gen]), accmean(&AVE_VP[gen]),
		accmean(&AVE_VA[gen]), accmean(&AVE_VD[gen]), accmean(&AVE_ID[gen]));

	fprintf(fgen, "\n*********** ID from homozygotes at gen 0 and the last generation ***********\n\n");

	if (neutral == 0)
	{
		fprintf(fgen, "gen 0  mean_FITN_0=%6.4f    mean_FITN_F_0=%6.4f    ID_F1_0=%6.4f\n", accmean(&AVE_mean_FITN_0), accmean(&AVE_mean_FITN_F_0), accmean(&AVE_ID_FITN_0)); 
		fprintf(fgen, "gen %d  mean_FITN=%6.4f    mean_FITN_F=%6.4f    ID_F1=%6.4f\n", generations, accmean(&AVE_mean_FITN), accmean(&AVE_mean_FITN_F), accmean(&AVE_ID_FITN)); 
	}
	else
	{
		fprintf(fgen, "gen 0  mean_PHEN_0=%6.4f    mean_PHEN_F_0=%6.4f    ID_F1_0=%6.4f\n", accmean(&AVE_mean_PHEN_0), accmean(&AVE_mean_PHEN_F_0), accmean(&AVE_ID_PHEN_0)); 
		fprintf(fgen, "gen %d  mean_PHEN=%6.4f    mean_PHEN_F=%6.4f    ID_F1=%6.4f\n", generations, accmean(&AVE_mean_PHEN), accmean(&AVE_mean_PHEN_F), accmean(&AVE_ID_PHEN)); 
	}
	fprintf(fgen, "\n*********** Frequencies at genBP ***********\n\n");

	fprintf(fgen, "FREQ 0.0   0.0-0.1   0.1-0.2   0.2-0.3   0.3-0.4   0.4-0.5   0.5-0.6   0.6-0.7   0.7-0.8   0.8-0.9   0.9-1.0   10\n");
	fprintf(fgen, "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
		accmean(&AVE_qi_0), accmean(&AVE_qi_0_1), accmean(&AVE_qi_1_2), accmean(&AVE_qi_2_3), accmean(&AVE_qi_3_4), accmean(&AVE_qi_4_5), accmean(&AVE_qi_5_6), accmean(&AVE_qi_6_7), accmean(&AVE_qi_7_8), accmean(&AVE_qi_8_9), accmean(&AVE_qi_9_10), accmean(&AVE_qi_10) );

	fprintf(fgen, "\n*********** Frequencies current generation ***********\n\n");

	fprintf(fgen, "MEAN F\ngen  LOCI   EHom    EHet   MAF     D2      r2      Fped    Fibd    Fvr1     Fvr2     Fyang1    Fyang2    FLH1    FLH2     Fhom\n");
	for(gen=0; gen<=generations; gen++)
	if (gen >= genBP)
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		gen, numLOCI[gen], accmean(&AVE_EHom[gen]), accmean(&AVE_EHet[gen]), accmean(&AVE_MAF[gen]), accmean(&AVE_D2[gen]), accmean(&AVE_r2[gen]), accmean(&AVE_Fped[gen]), accmean(&AVE_Fibd[gen]), accmean(&AVE_Fvr1[gen]), accmean(&AVE_Fvr2[gen]), accmean(&AVE_Fyang1[gen]), accmean(&AVE_Fyang2[gen]), accmean(&AVE_FLH1[gen]), accmean(&AVE_FLH2[gen]), accmean(&AVE_Fhom[gen]));

//	GRAPH
	fprintf(fGRAPH, " 00000  00000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_Fvr1[generations]), accmean(&AVE_Fvr2[generations]), accmean(&AVE_Fyang1[generations]), accmean(&AVE_Fyang2[generations]), accmean(&AVE_FLH1[generations]), accmean(&AVE_FLH2[generations]), accmean(&AVE_Fhom[generations]));

	fprintf(fgen, "SD F\ngen  LOCI   Fped    Fibd    Fvr1    Fvr2    Fyang1   Fyang2   FLH1    FLH2    Fhom\n");
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		generations, numLOCI[generations], sqrt(variance(&AVE_Fped[generations])), sqrt(variance(&AVE_Fibd[generations])), sqrt(variance(&AVE_Fvr1[generations])), sqrt(variance(&AVE_Fvr2[generations])), sqrt(variance(&AVE_Fyang1[generations])), sqrt(variance(&AVE_Fyang2[generations])), sqrt(variance(&AVE_FLH1[generations])), sqrt(variance(&AVE_FLH2[generations])), sqrt(variance(&AVE_Fhom[generations])));

	fprintf(fgen, "VAR F\ngen  LOCI   Fped    Fibd    Fvr1    Fvr2    Fyang1   Fyang2   FLH1    FLH2    Fhom\n");
	for(gen=0; gen<=generations; gen++)
	if (gen >= genBP)
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		generations, numLOCI[gen], accmean(&AVE_VFped[gen]), accmean(&AVE_VFibd[gen]), accmean(&AVE_VFvr1[gen]), accmean(&AVE_VFvr2[gen]), accmean(&AVE_VFyang1[gen]), accmean(&AVE_VFyang2[gen]), accmean(&AVE_VFLH1[gen]), accmean(&AVE_VFLH2[gen]), accmean(&AVE_VFhom[gen]));

//	GRAPH
	fprintf(fGRAPH, " 00000  00000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_VFvr1[generations]), accmean(&AVE_VFvr2[generations]), accmean(&AVE_VFyang1[generations]), accmean(&AVE_VFyang2[generations]), accmean(&AVE_VFLH1[generations]), accmean(&AVE_VFLH2[generations]), accmean(&AVE_VFhom[generations]));

	fprintf(fgen, "SD VAR F\ngen  LOCI   Fped    Fibd    Fvr1    Fvr2    Fyang1   Fyang2   FLH1    FLH2    Fhom\n");
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		generations, numLOCI[generations], sqrt(variance(&AVE_VFped[generations])), sqrt(variance(&AVE_VFibd[generations])), sqrt(variance(&AVE_VFvr1[generations])), sqrt(variance(&AVE_VFvr2[generations])), sqrt(variance(&AVE_VFyang1[generations])), sqrt(variance(&AVE_VFyang2[generations])), sqrt(variance(&AVE_VFLH1[generations])), sqrt(variance(&AVE_VFLH2[generations])));

	fprintf(fgen, "\nMEAN ID at generation %d with current frequencies \n", generations); 
	fprintf(fgen, "       ID     IDF1W    IDF1Q    Fvr1     Fvr2     Fyang1    Fyang2    FLH1     FLH2     Fhom\n"); 
	fprintf(fgen, "%3d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", generations, -accmean(&AVE_ID[generations]), -accmean(&AVE_ID_FITN), accmean(&AVE_ID_PHEN), accmean(&AVE_IDFvr1), accmean(&AVE_IDFvr2), accmean(&AVE_IDFyang1), accmean(&AVE_IDFyang2), accmean(&AVE_IDFLH1), accmean(&AVE_IDFLH2), accmean(&AVE_IDFhom));
	fprintf(fgen, "\nSD ID\n       ID    IDF1W    IDF1Q   Fvr1    Fvr2    Fyang1   Fyang2   FLH1   FLH2     Fhom\n"); 
	fprintf(fgen, "%3d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", generations, sqrt(variance(&AVE_ID[generations])), sqrt(variance(&AVE_ID_FITN)), sqrt(variance(&AVE_ID_PHEN)), sqrt(variance(&AVE_IDFvr1)), sqrt(variance(&AVE_IDFvr2)), sqrt(variance(&AVE_IDFyang1)), sqrt(variance(&AVE_IDFyang2)), sqrt(variance(&AVE_IDFLH1)), sqrt(variance(&AVE_IDFLH2)), sqrt(variance(&AVE_IDFhom)));

//	GRAPH
	fprintf(fGRAPH, "00000  00000  00000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_IDFvr1), accmean(&AVE_IDFvr2), accmean(&AVE_IDFyang1), accmean(&AVE_IDFyang2), accmean(&AVE_IDFLH1), accmean(&AVE_IDFLH2), accmean(&AVE_IDFhom));

	fprintf(fgen, "\n       Fibd    Fvr1    Fvr2    Fyang1   Fyang2   FLH1    FLH2    Fhom\n"); 
	fprintf(fgen, "Fped  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FpedFibd), accmean(&AVE_FpedFvr1), accmean(&AVE_FpedFvr2), accmean(&AVE_FpedFyang1), accmean(&AVE_FpedFyang2), accmean(&AVE_FpedFLH1), accmean(&AVE_FpedFLH2), accmean(&AVE_FpedFhom)); 
	fprintf(fgen, "Fibd          %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FibdFvr1), accmean(&AVE_FibdFvr2), accmean(&AVE_FibdFyang1), accmean(&AVE_FibdFyang2), accmean(&AVE_FibdFLH1), accmean(&AVE_FibdFLH2), accmean(&AVE_FibdFhom)); 
	fprintf(fgen, "Fvr1                  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr1Fvr2), accmean(&AVE_Fvr1Fyang1), accmean(&AVE_Fvr1Fyang2), accmean(&AVE_Fvr1FLH1), accmean(&AVE_Fvr1FLH2), accmean(&AVE_Fvr1Fhom)); 
	fprintf(fgen, "Fvr2                          %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr2Fyang1), accmean(&AVE_Fvr2Fyang2), accmean(&AVE_Fvr2FLH1), accmean(&AVE_Fvr2FLH2), accmean(&AVE_Fvr2Fhom)); 
	fprintf(fgen, "Fyang1                                %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fyang1Fyang2), accmean(&AVE_Fyang1FLH1), accmean(&AVE_Fyang1FLH2), accmean(&AVE_Fyang1Fhom)); 
	fprintf(fgen, "Fyang2                                        %7.4f %7.4f %7.4f\n", accmean(&AVE_Fyang2FLH1), accmean(&AVE_Fyang2FLH2), accmean(&AVE_Fyang2Fhom)); 
	fprintf(fgen, "FLH1                                                  %7.4f %7.4f\n", accmean(&AVE_FLH1FLH2), accmean(&AVE_FLH1Fhom)); 
	fprintf(fgen, "FLH2                                                          %7.4f\n\n", accmean(&AVE_FLH2Fhom)); 

//	GRAPH
	fprintf(fGRAPH, "  00000  00000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n\n",
		accmean(&AVE_FibdFvr1), accmean(&AVE_FibdFvr2), accmean(&AVE_FibdFyang1), accmean(&AVE_FibdFyang2), accmean(&AVE_FibdFLH1), accmean(&AVE_FibdFLH2), accmean(&AVE_FibdFhom));

	fprintf(fgen, "\nNeFped  NeFibd  NeFvr1 NeFvr2  NeFyang1 NeFyang2 NeFLH1 NeFLH1\n");
	fprintf(fgen, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_NeFped), accmean(&AVE_NeFibd), accmean(&AVE_NeFvr1), accmean(&AVE_NeFvr2), accmean(&AVE_NeFyang1), accmean(&AVE_NeFyang2), accmean(&AVE_NeFLH1), accmean(&AVE_NeFLH2), accmean(&AVE_NeFhom) );

	fprintf(fgen, "\n*********** Frequencies BP=%d generation ***********\n", genBP);

	fprintf(fgen, "\nMEAN F(BP)\ngen  LOCI   EHom    EHet   MAF     D2      r2      Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP\n");
	for(gen=0; gen<=generations; gen++)
	if (gen >= genBP)
	{
		fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			gen, numLOCI[gen], accmean(&AVE_EHomBP[gen]), accmean(&AVE_EHetBP[gen]), accmean(&AVE_MAF[gen]), accmean(&AVE_D2[gen]), accmean(&AVE_r2[gen]), accmean(&AVE_Fped[gen]), accmean(&AVE_Fibd[gen]), accmean(&AVE_Fvr1BP[gen]), accmean(&AVE_Fvr2BP[gen]), accmean(&AVE_Fyang1BP[gen]), accmean(&AVE_Fyang2BP[gen]), accmean(&AVE_FLH1BP[gen]), accmean(&AVE_FLH2BP[gen]));

		if ((gen == generations)&&(rep==replicates)) fprintf(fOUTM, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			gen, numLOCI[gen], accmean(&AVE_MAF[gen]), accmean(&AVE_D2[gen]), accmean(&AVE_r2[gen]), accmean(&AVE_Fped[gen]), accmean(&AVE_Fibd[gen]), accmean(&AVE_Fvr1BP[gen]), accmean(&AVE_Fvr2BP[gen]), accmean(&AVE_Fyang1BP[gen]), accmean(&AVE_Fyang2BP[gen]), accmean(&AVE_FLH1BP[gen]), accmean(&AVE_FLH2BP[gen]));
	}

//	GRAPH
	fprintf(fGRAPH, " %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_Fibd[generations]), accmean(&AVE_Fped[generations]), accmean(&AVE_Fvr1BP[generations]), accmean(&AVE_Fvr2BP[generations]), accmean(&AVE_Fyang1BP[generations]), accmean(&AVE_Fyang2BP[generations]), accmean(&AVE_FLH1BP[generations]), accmean(&AVE_FLH2BP[generations]));

	fprintf(fgen, "SD F(BP)\ngen  LOCI   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP\n");
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		generations, numLOCI[generations], sqrt(variance(&AVE_Fped[generations])), sqrt(variance(&AVE_Fibd[generations])), sqrt(variance(&AVE_Fvr1BP[generations])), sqrt(variance(&AVE_Fvr2BP[generations])), sqrt(variance(&AVE_Fyang1BP[generations])), sqrt(variance(&AVE_Fyang2BP[generations])), sqrt(variance(&AVE_FLH1BP[generations])), sqrt(variance(&AVE_FLH2BP[generations])));

	fprintf(fgen, "VAR F(BP)\ngen  LOCI   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP\n");
	for(gen=0; gen<=generations; gen++)
	if (gen >= genBP)
	{
		fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			gen, numLOCI[gen], accmean(&AVE_VFped[gen]), accmean(&AVE_VFibd[gen]), accmean(&AVE_VFvr1BP[gen]), accmean(&AVE_VFvr2BP[gen]), accmean(&AVE_VFyang1BP[gen]), accmean(&AVE_VFyang2BP[gen]), accmean(&AVE_VFLH1BP[gen]), accmean(&AVE_VFLH2BP[gen]));

		if ((gen == generations)&&(rep==replicates))	fprintf(fOUTV, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			gen, numLOCI[gen], accmean(&AVE_VFped[gen]), accmean(&AVE_VFibd[gen]), accmean(&AVE_VFvr1BP[gen]), accmean(&AVE_VFvr2BP[gen]), accmean(&AVE_VFyang1BP[gen]), accmean(&AVE_VFyang2BP[gen]), accmean(&AVE_VFLH1BP[gen]), accmean(&AVE_VFLH2BP[gen]));
	}

//	GRAPH
	fprintf(fGRAPH, " %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_VFibd[generations]), accmean(&AVE_VFped[generations]), accmean(&AVE_VFvr1BP[generations]), accmean(&AVE_VFvr2BP[generations]), accmean(&AVE_VFyang1BP[generations]), accmean(&AVE_VFyang2BP[generations]), accmean(&AVE_VFLH1BP[generations]), accmean(&AVE_VFLH2BP[generations]));

	fprintf(fgen, "SD VAR F(BP)\ngen  LOCI   Fped    Fibd    Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP\n");
	fprintf(fgen, "%3d  %d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		generations, numLOCI[generations], sqrt(variance(&AVE_VFped[generations])), sqrt(variance(&AVE_VFibd[generations])), sqrt(variance(&AVE_VFvr1BP[generations])), sqrt(variance(&AVE_VFvr2BP[generations])), sqrt(variance(&AVE_VFyang1BP[generations])), sqrt(variance(&AVE_VFyang2BP[generations])), sqrt(variance(&AVE_VFLH1BP[generations])), sqrt(variance(&AVE_VFLH2BP[generations])));

	fprintf(fgen, "\nMEAN ID(BP) at generation %d with frequencies at generation %d\n", generations, genBP); 
	fprintf(fgen, "       ID      IDF1W    IDF1Q    Fped    Fibd    Fvr1BP   Fvr2BP   Fyang1BP  Fyang2BP  FLH1BP   FLH2BP\n"); 
	fprintf(fgen, "%3d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", generations, -accmean(&AVE_ID[genBP]), accmean(&AVE_ID_FITN_0), -accmean(&AVE_ID_PHEN_0), accmean(&AVE_IDFped), accmean(&AVE_IDFibd), accmean(&AVE_IDFvr1BP), accmean(&AVE_IDFvr2BP), accmean(&AVE_IDFyang1BP), accmean(&AVE_IDFyang2BP), accmean(&AVE_IDFLH1BP), accmean(&AVE_IDFLH2BP));

	fprintf(fgen, "\nSD ID(BP)\n      ID    IDF1W    IDF1Q    Fped   Fibd   Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP  FLH2BP\n"); 
	fprintf(fgen, "%3d  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n", generations, sqrt(variance(&AVE_ID[generations])), sqrt(variance(&AVE_ID_FITN_0)), sqrt(variance(&AVE_ID_PHEN_0)), sqrt(variance(&AVE_IDFped)), sqrt(variance(&AVE_IDFibd)), sqrt(variance(&AVE_IDFvr1BP)), sqrt(variance(&AVE_IDFvr2BP)), sqrt(variance(&AVE_IDFyang1BP)), sqrt(variance(&AVE_IDFyang2BP)), sqrt(variance(&AVE_IDFLH1BP)), sqrt(variance(&AVE_IDFLH2BP)));

//	GRAPH
	fprintf(fGRAPH, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		-accmean(&AVE_ID[generations]), accmean(&AVE_IDFibd), accmean(&AVE_IDFped), accmean(&AVE_IDFvr1BP), accmean(&AVE_IDFvr2BP), accmean(&AVE_IDFyang1BP), accmean(&AVE_IDFyang2BP), accmean(&AVE_IDFLH1BP), accmean(&AVE_IDFLH2BP));

	fprintf(fgen, "\nCorrelations at generation %d with frequencies of generation %d\n", generations, genBP); 

	fprintf(fgen, "\n        FibdBP  Fvr1BP  Fvr2BP  Fyang1BP Fyang2BP FLH1BP FLH2BP\n"); 
	fprintf(fgen, "Fped  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FpedFibdBP), accmean(&AVE_FpedFvr1BP), accmean(&AVE_FpedFvr2BP), accmean(&AVE_FpedFyang1BP), accmean(&AVE_FpedFyang2BP), accmean(&AVE_FpedFLH1BP), accmean(&AVE_FpedFLH2BP)); 
	fprintf(fgen, "Fibd          %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_FibdFvr1BP), accmean(&AVE_FibdFvr2BP), accmean(&AVE_FibdFyang1BP), accmean(&AVE_FibdFyang2BP), accmean(&AVE_FibdFLH1BP), accmean(&AVE_FibdFLH2BP)); 
	fprintf(fgen, "Fvr1BP                %7.4f %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr1Fvr2BP), accmean(&AVE_Fvr1Fyang1BP), accmean(&AVE_Fvr1Fyang2BP), accmean(&AVE_Fvr1FLH1BP), accmean(&AVE_Fvr1FLH2BP)); 
	fprintf(fgen, "Fvr2BP                        %7.4f %7.4f %7.4f %7.4f\n", accmean(&AVE_Fvr2Fyang1BP), accmean(&AVE_Fvr2Fyang2BP), accmean(&AVE_Fvr2FLH1BP), accmean(&AVE_Fvr2FLH2BP)); 
	fprintf(fgen, "Fyang1BP                              %7.4f %7.4f %7.4f\n", accmean(&AVE_Fyang1Fyang2BP), accmean(&AVE_Fyang1FLH1BP), accmean(&AVE_Fyang1FLH2BP)); 
	fprintf(fgen, "Fyang2BP                               	      %7.4f %7.4f\n", accmean(&AVE_Fyang2FLH1BP), accmean(&AVE_Fyang2FLH2BP)); 
	fprintf(fgen, "FLH1BP                                                %7.4f\n", accmean(&AVE_FLH1FLH2BP)); 

//	GRAPH
	fprintf(fGRAPH, " 00000  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_FpedFibdBP), accmean(&AVE_FibdFvr1BP), accmean(&AVE_FibdFvr2BP), accmean(&AVE_FibdFyang1BP), accmean(&AVE_FibdFyang2BP), accmean(&AVE_FibdFLH1BP), accmean(&AVE_FibdFLH2BP));

	fprintf(fgen, "\nNeFpedBP  NeFibdBP  NeFvr1BP NeFvr2BP  NeFyang1BP NeFyang2BP NeFLH1BP NeFLH1BP\n");
	fprintf(fgen, "%6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
		accmean(&AVE_NeFpedBP), accmean(&AVE_NeFibdBP), accmean(&AVE_NeFvr1BP), accmean(&AVE_NeFvr2BP), accmean(&AVE_NeFyang1BP), accmean(&AVE_NeFyang2BP), accmean(&AVE_NeFLH1BP), accmean(&AVE_NeFLH2BP) );

	fclose(fGRAPH);
	fclose(fgen);
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

