// TEST_PROB3.c (14/08/2022)

#include "libhdr"
#define NN 2001
#define PP 10001

int i, j, k, rnd, REPS, PERMUTS;
double w, TEST[NN], ID[NN], F1[NN], ibd[NN], ped[NN], vr1[NN], vr2[NN];
double yan1[NN], yan2[NN], lh1[NN], lh2[NN], hom[NN], R100[NN], R1000[NN], R5000[NN];
double DIFF[PP][15], Prob_menor_0_[15], Prob_mayor_0_[15];

struct acc OBS_DIFF[NN], AVE_DIFF[15];

FILE *fdat, *fout;

main()
{
	fout = fopen ("outfile","w");

	readfiles();
	differences();
	permutations();
	probabilities();
	printing();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	tracestart();
	getseed();
	getintandskip("Number of replicates :",&REPS,1,2000);
	getintandskip("Number of bootstraps :",&PERMUTS,1,10000);

	fdat = fopen ("data","r");
	for (i=1; i<=REPS; i++)
	{
		fscanf(fdat,"%lf", &w);
		ID[i] = w;
		fscanf(fdat,"%lf", &w);
		F1[i] = w;
		fscanf(fdat,"%lf", &w);
		ibd[i] = w;
		fscanf(fdat,"%lf", &w);
		ped[i] = w;
		fscanf(fdat,"%lf", &w);
		vr1[i] = w;
		fscanf(fdat,"%lf", &w);
		vr2[i] = w;
		fscanf(fdat,"%lf", &w);
		yan1[i] = w;
		fscanf(fdat,"%lf", &w);
		yan2[i] = w;
		fscanf(fdat,"%lf", &w);
		lh1[i] = w;
		fscanf(fdat,"%lf", &w);
		lh2[i] = w;
		fscanf(fdat,"%lf", &w);
		hom[i] = w;
		fscanf(fdat,"%lf", &w);
		R100[i] = w;
		fscanf(fdat,"%lf", &w);
		R1000[i] = w;
		fscanf(fdat,"%lf", &w);
		R5000[i] = w;
	}
	fclose(fdat);

	for (i=1; i<=REPS; i++)	TEST[i] = ibd[i];

/*	for (i=1; i<=3; i++)
	{
		fprintf(fout, "ID=%f F1=%f ibd=%f ped=%f vr1=%f vr2=%f yan1=%f yan2=%f lh1=%f lh2=%f hom=%f R100=%f R1000=%f R5000=%f\n",
			ID[i], F1[i], ibd[i], ped[i], vr1[i], vr2[i], yan1[i], yan2[i], lh1[i], lh2[i], hom[i], R100[i], R1000[i], R5000[i]);
	}
*/

	if (tracelevel != 0)
	{
		fprintf(fout, "\nINITIAL\n");
		for (i=1; i<=REPS; i++)	fprintf(fout, "ibd=%f yan2=%f\n", ibd[i], yan2[i]);
		fprintf(fout, "\n");
	}

	writeseed();
	return(0);
}

/* **************************************************************************** */

differences()
{
	for (i=1; i<=REPS; i++)
	{
		accum(&OBS_DIFF[1], (TEST[i] - ID[i]));
		accum(&OBS_DIFF[2], (TEST[i] - F1[i]));
		accum(&OBS_DIFF[3], (TEST[i] - ibd[i]));
		accum(&OBS_DIFF[4], (TEST[i] - ped[i]));
		accum(&OBS_DIFF[5], (TEST[i] - vr1[i]));
		accum(&OBS_DIFF[6], (TEST[i] - vr2[i]));
		accum(&OBS_DIFF[7], (TEST[i] - yan1[i]));
		accum(&OBS_DIFF[8], (TEST[i] - yan2[i]));
		accum(&OBS_DIFF[9], (TEST[i] - lh1[i]));
		accum(&OBS_DIFF[10], (TEST[i] - lh2[i]));
		accum(&OBS_DIFF[11], (TEST[i] - hom[i]));
		accum(&OBS_DIFF[12], (TEST[i] - R100[i]));
		accum(&OBS_DIFF[13], (TEST[i] - R1000[i]));
		accum(&OBS_DIFF[14], (TEST[i] - R5000[i]));
	}

	if (tracelevel != 0)
	{
		fprintf(fout, "DIFFERENCE\n");
		fprintf(fout, "DIFF_ibd_yan2=%f\n", accmean(&OBS_DIFF[8]));
	}

	return(0);
}

/* **************************************************************************** */

permutations()
{
	fprintf(fout, "\nPERMUTS\n");

	for (i=1; i<=PERMUTS; i++)
	{
		if (tracelevel != 0)	for (j=1; j<=PERMUTS; j++)	fprintf(fout, "i=%d TEST(%d)=%f yan2(%d)=%f\n", i, j, TEST[j], j, yan2[j]);

		for (j=1; j<=PERMUTS; j++)
		{
			rnd = (int)(uniform() * REPS) + 1;

			DIFF[i][1] += (TEST[rnd] - ID[rnd]) / PERMUTS;
			DIFF[i][2] += (TEST[rnd] - F1[rnd]) / PERMUTS;
			DIFF[i][3] += (TEST[rnd] - ibd[rnd]) / PERMUTS;
			DIFF[i][4] += (TEST[rnd] - ped[rnd]) / PERMUTS;
			DIFF[i][5] += (TEST[rnd] - vr1[rnd]) / PERMUTS;
			DIFF[i][6] += (TEST[rnd] - vr2[rnd]) / PERMUTS;
			DIFF[i][7] += (TEST[rnd] - yan1[rnd]) / PERMUTS;
			DIFF[i][8] += (TEST[rnd] - yan2[rnd]) / PERMUTS;
			DIFF[i][9] += (TEST[rnd] - lh1[rnd]) / PERMUTS;
			DIFF[i][10] += (TEST[rnd] - lh2[rnd]) / PERMUTS;
			DIFF[i][11] += (TEST[rnd] - hom[rnd]) / PERMUTS;
			DIFF[i][12] += (TEST[rnd] - R100[rnd]) / PERMUTS;
			DIFF[i][13] += (TEST[rnd] - R1000[rnd]) / PERMUTS;
			DIFF[i][14] += (TEST[rnd] - R5000[rnd]) / PERMUTS;
		}

		if (tracelevel != 0)
		{
			fprintf(fout, "\nDIFF PERMUTS\n");
			fprintf(fout, "i=%d DIFF_ibd_yan2=%f\n", i, DIFF[i][8]);
			fprintf(fout, "\n\n");
		}
	}

	// 95% intervals

	for (k=1; k<=14; k++)
	for (i=1; i<PERMUTS; i++)
	for (j=i+1; j<=PERMUTS; j++)
	{
		if (DIFF[j][k] < DIFF[i][k])
		{
			w = DIFF[i][k]; DIFF[i][k] = DIFF[j][k]; DIFF[j][k] = w;
		}
	}

	if (tracelevel != 0)
	{
		for (i=1; i<=PERMUTS; i++)
		fprintf(fout, "DIFF[%d][8]=%f\n", i, DIFF[i][8]);
	}

	for (i=1; i<=PERMUTS; i++)
	for (k=1; k<=14; k++)	accum(&AVE_DIFF[k], DIFF[i][k]);

	if (tracelevel != 0)
	{
		fprintf(fout, "\nAVE DIFF PERMUTS\n");
		fprintf(fout, "AVE_DIFF_ibd_yan2=%f   lower=%f   upper=%f\n", accmean(&AVE_DIFF[8]), DIFF[(int)(PERMUTS*0.025)+1][8], DIFF[(int)(PERMUTS*0.975)+1][8]);
		fprintf(fout, "\n\n");
	}
}

/* **************************************************************************** */

probabilities()
{
	for (k=1; k<=14; k++)
	{
		if (tracelevel != 0)
		{
			fprintf(fout, "OBS_DIFF[%d]=%f   AVE_DIFF=%f   lower=%f   upper=%f\n", k, accmean(&OBS_DIFF[k]), accmean(&AVE_DIFF[k]), DIFF[(int)(PERMUTS*0.025)+1][k], DIFF[(int)(PERMUTS*0.975)+1][k]);
		}

		if ((accmean(&OBS_DIFF[k]) - DIFF[i][k]) < 0.0)	Prob_menor_0_[k] ++;
		else										Prob_mayor_0_[k] ++;
	}

	for (k=1; k<=14; k++)
	for (i=1; i<PERMUTS; i++)
	{
		if (DIFF[i][k] < 0.0)	Prob_menor_0_[k] ++;
		else					Prob_mayor_0_[k] ++;
	}

	return(0);
}

/* **************************************************************************** */

printing()
{
	for (k=1; k<=14; k++)
	fprintf(fout, "k=%d  OBS_DIFF=%f   AVE_DIFF=%f   lower=%f   upper=%f\n", k, accmean(&OBS_DIFF[k]), accmean(&AVE_DIFF[k]), DIFF[(int)(PERMUTS*0.025)+1][k], DIFF[(int)(PERMUTS*0.975+1)][k]);

	for (k=1; k<=14; k++)
	fprintf(fout,"%d    %f   %f\n", k, Prob_menor_0_[k]/PERMUTS, Prob_mayor_0_[k]/PERMUTS);

	return(0);
}

/* **************************************************************************** */

