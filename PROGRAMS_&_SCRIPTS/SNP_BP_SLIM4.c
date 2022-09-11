// SNP_BP_SLIM4.c

#include "libhdr"

#define NN 300000 // Maximum 300000 SNPs segregating
#define CC 1001   // Maximum N=1000 individuals

int i, j, s, ss, m, x, b, nind, numSNP;
unsigned long long int pos[NN], xx, a;
int mut[NN], crom[NN], num_mut_crom[CC], snp[CC][NN];
double w, ps[NN], ef[NN], h[NN], q[NN];
char ch;

FILE *fdat, *fped, *fmap, *fallsnps;

main()
{
	getseed();
	getintandskip("NIND :",&nind,10,1000);

	readfiles();
	PLINK_files();
	writeseed();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	fped = fopen ("dataBP.ped","w");
	fmap = fopen ("dataBP.map","w");
	fallsnps = fopen ("list_allsnps","w");

	// ********** Read slimout to get SNP positions, effects and frequencies ********** 

	fdat = fopen ("slimout","r");

	// ********** gets the position, effects and frequencies of all mutations simulated (numSNP) **********

	lookfortext("Mutations:");

	while (!feof(fdat))
	{
		s ++;
		fscanf(fdat,"%d", &x);
		mut[s] = x;
		fscanf(fdat,"%d", &x);
		for (j=1; j<=4; j++)	fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%llu", &xx);
		pos[s] = xx;
		if (pos[s] == pos[s-1])
		{
			pos[s] = pos[s-1] + 1;
		}
		fscanf(fdat,"%lf", &w);
		ps[s] = w;
		ef[s] = fabs(w);
		fscanf(fdat,"%lf", &w);
		h[s] = w;
		if (ps[s] == 0.0)	h[s] = 0.0;
		for (j=1; j<=4; j++)	fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);
		fscanf(fdat,"%d", &x);
		q[s] = x / (2.0*nind);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		if (ch != 'I')	ungetc(ch, fdat);
		else		break;
	}
	numSNP = s;

	fclose(fdat);

	// ********** reorder by genome position **********

	for (s=1; s<numSNP; s++)
	for (ss=s; ss<=numSNP; ss++)
	{
		if (pos[ss] < pos[s])
		{
			a=pos[s]; pos[s]=pos[ss]; pos[ss]=a;
			b=mut[s]; mut[s]=mut[ss]; mut[ss]=b;
			w=ps[s]; ps[s]=ps[ss]; ps[ss]=w;
			w=ef[s]; ef[s]=ef[ss]; ef[ss]=w;
			w=h[s]; h[s]=h[ss]; h[ss]=w;
			w=q[s]; q[s]=q[ss]; q[ss]=w;
		}
	}

	// ********** gets mutations **********

	fdat = fopen ("slimout","r");
	lookfortext("Genomes:");

	fscanf(fdat,"%c", &ch);

	for (i=1; i<=(2*nind);i++)
	{
		m = 0;

		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%d", &x);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		fscanf(fdat,"%c", &ch);
		if (ch == '\n')	goto next;
		else	ungetc(ch, fdat);
		while (!feof(fdat))
		{
			fscanf(fdat,"%c", &ch);

			if (ch == '\n')	break;
			else
			{
				ungetc(ch, fdat);
				m ++;
				fscanf(fdat,"%d", &x);
				crom[m] = x;
			}			
		}
		num_mut_crom[i] = m;

		// ********** Identify SNPs **********

		for (s=1; s<=numSNP; s++)
		{
			snp[i][s] = 0;

			for (m=1; m<=num_mut_crom[i]; m++)
			{
				if (crom[m] == mut[s])
				{
					snp[i][s] = 1;
					break;
				}
			}
		}

		next: /***/;
	}

	fclose(fdat);

	return(0);
}

/* **************************************************************************** */

PLINK_files()
{
	// ********** PLINK data.map file **********

	for (s=1; s<=numSNP; s++)	fprintf(fmap,"%llu SNP%d 0 %llu\n", (pos[s]/125000000)+1, s, pos[s]);

	// ********** PLINK data.ped file **********

	for (i=1; i<(2*nind);i+=2)
	{
		fprintf(fped,"1 IND%d 0 0 1 0   ", (i+1)/2);

		for (s=1; s<=numSNP; s++)
		{
			if ((snp[i][s]==0) && (snp[i+1][s]==0))			fprintf(fped,"1 1  ");
			else if ((snp[i][s]==0) && (snp[i+1][s]==1))		fprintf(fped,"1 2  ");
			else if ((snp[i][s]==1) && (snp[i+1][s]==0))		fprintf(fped,"2 1  ");
			else if ((snp[i][s]==1) && (snp[i+1][s]==1))		fprintf(fped,"2 2  ");
		}
	}
	fprintf(fped,"\n");

	// List all SNPs

	fprintf(fallsnps,"%d\n", numSNP);
	for (s=1; s<=numSNP; s++)	fprintf(fallsnps,"%llu  %llu  %f  %f  %f  %f\n", (pos[s]/125000000)+1, pos[s], ps[s], ef[s], h[s], q[s]);

	return(0);
}

/* **************************************************************************** */

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
      c = getc(fdat);

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

