// SNP_BP_SLIM3.c

#include "libhdr"

#define NN 200000 // Maximum 200000 SNPs segregating
#define CC 1001   // Maximum N=1000 individuals
#define MC 120000 // Maximum 120000 SNPs per chromosome

int i, j, s, ss, m, x, b, nind, HOM1, HOM2, crossover;
unsigned long long int pos[NN], xx, a;
int mut[NN], crom[CC][MC], num_mut_crom[CC], chromosomeP[CC][NN][2];
double w, ps[NN], ef[NN], h[NN], q[NN], FITN[CC];
double FITNparents[CC];
int numSNPchip, SNPchip[NN], numSNP;
char ch;

FILE *fdat, *fgen, *fout, *fPpos, *fPphen, *fallsnps, *fLISTQTL;

main()
{
	getseed();
	getintandskip("NIND :",&nind,10,1000);

	readfiles();
	PLINK_files();
	printing();
	writeseed();

	return(0);
}

/* **************************************************************************** */

readfiles()
{
	fgen = fopen ("dataBP.ped","w");
	fout = fopen ("outfileSLIM","w");

	fPpos = fopen ("dataBP.map","w");
	fPphen = fopen ("qtBP.phe","w");

	fallsnps = fopen ("list_allsnps","w");
	fLISTQTL = fopen ("list_qtls","w");

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
				crom[i][m] = x;
			}			
		}
		num_mut_crom[i] = m;

		next: /***/;
	}

	fclose(fdat);

	return(0);
}

/* **************************************************************************** */

PLINK_files()
{
	double sum, sum2;

	for (i=1; i<(2*nind);i+=2)
	{
		if (i%2 != 0)
		{
			fprintf(fgen,"1 IND%d 0 0 1 0   ", (i+1)/2);
			FITN[(i+1)/2] = 1.0;
		}
		for (s=1; s<=numSNP; s++)
		{
			if (i == 1)
			{
				// PLINK POSFILE
				fprintf(fPpos,"%llu SNP%d 0 %llu\n", (pos[s]/125000000)+1, s, pos[s]);
			}
			HOM1 = 0; HOM2 = 0;

			for (m=1; m<=num_mut_crom[i]; m++)
			{
				if (crom[i][m] == mut[s])
				{
					HOM1 = 1;
					chromosomeP[(i+1)/2][s][1] = 1;
					break;
				}
			}
			for (m=1; m<=num_mut_crom[i+1]; m++)
			{
				if (crom[i+1][m] == mut[s])
				{
					HOM2 = 1;
					chromosomeP[(i+1)/2][s][2] = 1;
					break;
				}
			}

			if ((HOM1==0) && (HOM2==0))		fprintf(fgen,"1 1  ");
			else if ((HOM1==1) && (HOM2==0))
			{
				fprintf(fgen,"2 1  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);
			}
			else if ((HOM1==0) && (HOM2==1))
			{
				fprintf(fgen,"1 2  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]*h[s]);	
			}
			else
			{
				fprintf(fgen,"2 2  ");
				FITN[(i+1)/2] *= (1.0 + ps[s]);
			}
		}
		fprintf(fgen,"\n");

		if (i%2 != 0)
		{
			// PLINK PHENOFILE
			fprintf(fPphen,"1 IND%d %f\n", (i+1)/2, FITN[(i+1)/2]);
		}
	}

	// List all SNPs
//	fprintf(fallsnps,"chr  pos  s  a  h  q\n");
	fprintf(fallsnps,"%d\n", numSNP);
	for (s=1; s<=numSNP; s++)	fprintf(fallsnps,"%llu  %llu  %f  %f  %f  %f\n", (pos[s]/125000000)+1, pos[s], ps[s], ef[s], h[s], q[s]);

	//FILE WITH QTLs
	fprintf(fLISTQTL,"SNP  pos  s  a  h  q\n");
	for (s=1; s<=numSNP; s++)
	if (ps[s] != 0)
	fprintf(fLISTQTL,"%d  %llu  %f  %f  %f  %f\n", s, pos[s], ps[s], ef[s], h[s], q[s]);

	return(0);
}

/* **************************************************************************** */

printing()
{
	int hom;

	fprintf(fout,"Positions of SNPs\n\n");
	for (s=1; s<=numSNP; s++)
	{
		fprintf(fout,"SNP%d  mut=%d  pos=%llu  ps=%f  h=%f  q=%f\n", s, mut[s], pos[s], ps[s], h[s], q[s]);
	}
	fprintf(fout,"\n");

	for (i=1; i<=(2*nind);i++)
	{
		fprintf(fout,"CHROM %d  ", i);

		for (m=1; m<=num_mut_crom[i]; m++)
		fprintf(fout,"%d ", crom[i][m]);

		fprintf(fout,"\n");
	}

	fprintf(fout,"\n\Chromosome of parents\n\n");
	for (i=1; i<=nind;i++)
	for (hom=1; hom<=2;hom++)
	{
		fprintf(fout,"ind %d chromosomeP %d   ", i, hom);
		for (s=1; s<=numSNP; s++)
		{
			fprintf(fout,"%d ", chromosomeP[i][s][hom]);
		}
		fprintf(fout,"\n");
	}

	return(0);
}

/* ********************************************************************************************* */

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

