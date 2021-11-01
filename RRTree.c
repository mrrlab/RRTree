/******************************************************************************

Program written by Marc Robinson-Rechavi, implementing methods used and
developped in the article:

Robinson M, Gouy M, Gautier C, Mouchiroud D (1998) "Sensitivity of the relative
-rate test to taxonomic sampling", Mol Biol Evol 15: 1091-1098.

Many functions relating to distance computation or to manipulation of trees
were originally written by Nicolas Galtier for the software PHYLO_WIN:

Galtier N, Gouy M, Gautier C (1996) "SEAVIEW and PHYLO_WIN: two graphic tools
for sequence alignment and molecular phylogeny", Comp Appl Biosc 12: 543-548.

Many thanks to Nicolas for letting me use his code.

An application note describing RRTree has been published:

Robinson-Rechavi M and Huchon D (2000) "RRTree: Relative-Rate Tests between
groups of sequences on a phylogenetic tree", Bioinformatics 16, 296-297.

Marc Robinson-Rechavi - 31/july/2001

Laboratoire de Biologie Moleculaire et Cellulaire
Ecole Normale Superieure de Lyon
46, Allee d'Italie
69364 Lyon Cedex 07
FRANCE

tel    : +33 - 4 72 72 86 85
fax    : +33 - 4 72 72 80 80
e-mail : marc.robinson@ens-lyon.fr

Source codes are yours to use and modify freely as long as it's not for
commercial purposes, and that I'm credited.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define MAXLEN 2000
#define NAMELEN 100
#define NAMELEN_CLUSTAL 16
#define DNA_POISSON 3./4.
#define PROT_POISSON 19./20.
#define PRECISION 0.0000001	/* for computation of the exact probability */
#define PI_CST 3.1415926
#define NOT_IN_FILE -2

#define WINDOWS 0

struct options
	{
	int seqtype;
	int topology;
	int first;
	int code_mt;
	int compute_cds[5];
	int compute_nc;
	int support_limit;
	int several_sequences;
	int verbose;
	};

struct sequences
	{
	char *seq;
	char *name;
	char *com;
	int lg;
	};

struct parameters_lwl
	{
	double ks;
	double vks;
	double ka;
	double vka;
	double ls;
	double la;
	double a0;
	double a2;
	double a4;
	double b0;
	double b2;
	double b4;
	double vb4;
	double l0;
	double l2;
	double l4;
	double as;	/*synonymous transitions*/
	double vas;
	double ba;	/*non synonymous transversions*/
	double vba;
	};

struct parameters_k
	{
	double k;
	double vk;
	double l;
	double a;	/*if Kimura 2p, transitions*/
	double b;	/*if Kimura 2p, transversions*/
	};

enum sequence_types {NC, CDS, PROT};	/* Non Coding, CoDing Sequence, PROTein */
enum choice_compute_cds {ks, ka, as, ba, b4};
enum choice_compute_nc {jc, k2, tn6};
enum file_formats {mase, fasta, phylip, clustal, gde, nexus, lintre, phyltest, mega};

/**************************** FUNCTION PROTOTYPES ****************************/

void Treatement_cds(FILE *out_text, FILE *out_table, char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin, struct options *choices, char ** lineage_name);
void Treatement_nc(FILE *out_text, FILE *out_table, char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin, struct options *choices, char ** lineage_name);

struct parameters_lwl Lwl(struct sequences s[2], struct options *choices);
int Num(char *cod, struct options *choices);
void PreFastlwl(double **, double **, double **, double **, double **, double **, double **, double **, double **, double **, struct options *choices);
int Fastlwl(char **seq, int nbseq, int lgseq, double **tti0, double **tti1, double **tti2, double **ttv0, double **ttv1, double **ttv2, double **tl0, double **tl1, double **tl2, struct parameters_lwl *results, struct options *choices);

struct parameters_k Kim2(struct sequences *to_treat);
struct parameters_k Poisson(struct sequences *to_treat, double fraction);

struct parameters_k NC_dist(struct sequences *to_treat, int method);
double NC_cov(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34, int method);

double *GC(struct sequences s, struct options *choices);

double CovKs(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34);
double CovKa(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34);
double CovAs(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34);
double CovBa(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34);
double CovB4(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34);
double CovK2(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34);
double CovPoisson(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34, double fraction);

double Exact_probability(double obs_value);
double Divide(double a, double b);

double **Weighting(char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin);
void Polytomies(char *text, int **tab, char **names, int nb_otu);
int Unroot(int** treer, int** treenr, int notu);
void Root_bi(int** tree, int** r_tree, int bi, int notu);
void Root(int** tree, int** r_tree, int otu, int notu);
void Remove_column(int nb_otu, int **tab_tree);
void Remove_line(int remove, int nb_otu, int **tab_tree, char **names_tree);
int Tree_ctot(char *input, int **tree, double *lgbi, double *lgbp, double* bootstrap, char **name, int *root);
int Count_otu(char *tree);
int Is_text_zero(char *text, int start);
char *Unsignificant_branches_to_zero(char *tree_in, int limit);
char *Make_unsignificant_dichotomies(char *tree_in, int previous_nb_add);
char *Read_tree(FILE *ftree, struct options *choices);
void Ctree_root(char* ctree);
int Rooted(char* carbre);

FILE *Open_file(char mode[], char *name);
int Read_mase(FILE *f, int nb, struct sequences *s, int choice_com);
int Read_fasta_gde(FILE *f, int nb, struct sequences *s, char before_name);
int Read_phylip(FILE *f, int to_read, int read, struct sequences *s, int lg_name);
int Read_clustal(FILE *f, int to_read, int read, struct sequences *s);
int Read_nexus(FILE *f, int to_read, int read, struct sequences *s);
int Read_lintre(FILE *f, int nb, struct sequences *s);
int Read_phyltest(FILE *f, int read, int to_read, struct sequences *s);
int Read_mega(FILE *f, int to_read, int read, struct sequences *s);
int What_file_format(FILE *f, FILE *in_command_file);
void What_sequence_type(char *seq, struct options *choices, FILE *in_command_file);
int Get_full_line(char *line, int maxlen, FILE *f);
int Nbseq_mase(FILE *f);
int Nbseq_fasta_gde_phyltest(FILE *f, char before_name);
int Nbseq_phylip(FILE *f);
int Nbseq_clustal(FILE *f);
int Nbseq_nexus(FILE *f);
int Nbseq_lintre(FILE *f);
int Nbseq_mega(FILE *f);

int Read_in_command_file_number(FILE *f, char *search);
char *Read_in_command_file_string(FILE *f, char *search);

struct sequences Copy_sequence(struct sequences original);
void Free_sequences(struct sequences *s, int nb);

int Answer(int default_answer);
void Leave_program(char *message);
char *Read_line(char *line, int maxlen, FILE *f);

void Allgap(struct sequences **lineage, int *nb, int n_lineages, char *charlist);
struct sequences Nogap_cds(struct sequences in, char *charlist);
struct sequences Nogap(struct sequences in, char *charlist);

void *Allocation(int nb, size_t taille, char *here);
void *De_allocation(void *p);

/* a few Nicolas Galtier function prototypes */
int catsite(char c1, char c2, char c3, int i, struct options *choices);
char transf(char nt1, char nt2, struct options *choices);
void titv1(char *cod1, char *cod2, double weight, double *ti, double *tv, double* l, struct options *choices);
void titv2(char *cod1, char *cod2, double *ti, double *tv, double* l, int *aa, double **rl, int* pos, struct options *choices);
void titv3(char *cod1, char *cod2, double *ti, double *tv, double* l, int *aa, double **rl, struct options *choices);
int retder(int *list);
void aj(int *list, int nb);

/********************************** MAIN *************************************/

main(int argc, char *argv[])
{
FILE *in, *out_text, *out_table, *ftree, *in_command_file, *out_command_file;
struct sequences *read, *p;
struct sequences **lineage;	/*lineage[0] = outgroup*/
int *nb;
int the_end=0, status, file_format, lg_name_phylip;
int i, j, k, n;
char bases[]={'A', 'C', 'G', 'T', 'U', 0};
char aminoacids[]={'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 0};
int n_lineages, num_lin[3];
struct options *choices;
char *read_name, *name_file, line[MAXLEN], **lineage_name;
char *temp_string, *tree_parenthesis;

printf("\n\n                                  RRTree\n");
printf("\n                    a program for relative-rate tests\n");
printf("\n                                 version 1.1.11\n\n");
printf("Every time you have a choice to be answered by numbers, the choice between brackets [] is the default_answer, chosen if you simply push the return key.\n");
printf("If -1 is not a suggested answer, it will lead you to exit the program immediately.\n\n");

in_command_file=NULL;
out_command_file=NULL;
name_file=NULL;

switch (argc)
	{
	case 3:
		Leave_program("The 'default choices' option is no longer supported, please check about command files in RRTree.doc");
	case 2:
		if (argv[1][0]=='-')
			Leave_program("The 'default choices' option is no longer supported, please check about command files in RRTree.doc");
		else
			{
			name_file=Allocation(NAMELEN, sizeof(char), "name_file");
			strcpy(name_file, argv[1]);
			in_command_file=Open_file("r", name_file);
			name_file=De_allocation(name_file);
			break;
			}
	default:
		printf("command file [none] : ");
		in_command_file=Open_file("r", name_file);
		if (in_command_file==stdin) in_command_file=NULL;
	}

if (in_command_file==NULL)
	{
	printf("To avoid answering many questions next time, you may use a command file (see RRT_doc.txt).\n");
	printf("Do you want RRTree to automatically build a command file for you ([1]/0) ? ");
	if (Answer(1))
		{
		printf("command file name : ");
		out_command_file=Open_file("w", NULL);
		fprintf(out_command_file, "# automatically generated command file for RRTree\n\n");
		}
	}

name_file=Read_in_command_file_string(in_command_file, "aligned");
if (name_file==NULL)
	{
	name_file=Allocation(NAMELEN, sizeof(char), "name_file");
	read_name=Allocation(NAMELEN, sizeof(char), "read_name");

	printf("file of aligned sequences: ");
	Read_line(read_name, NAMELEN, stdin);
	sscanf(read_name, "%s", name_file);

	if (out_command_file) fprintf(out_command_file, "file of aligned sequences: %s\n", read_name);

	read_name=De_allocation(read_name);
	}

in=Open_file("r", name_file);

file_format=What_file_format(in, in_command_file);

if (out_command_file)
	switch (file_format)
		{
		case mase:
			fprintf(out_command_file, "format: mase\n");
			break;
		case fasta:
			fprintf(out_command_file, "format: fasta\n");
			break;
		case gde:
			fprintf(out_command_file, "format: gde\n");
			break;
		case phylip:
			fprintf(out_command_file, "format: phylip\n");
			break;
		case clustal:
			fprintf(out_command_file, "format: clustal\n");
			break;
		case nexus:
			fprintf(out_command_file, "format: nexus\n");
			break;
		case lintre:
			fprintf(out_command_file, "format: lintre\n");
			break;
		case phyltest:
			fprintf(out_command_file, "format: phyltest\n");
			break;
		case mega:
			fprintf(out_command_file, "format: mega\n");
			break;
		default:
			Leave_program("sorry, file format not implemented.");
		}

if (file_format==phylip)
	{
	lg_name_phylip=Read_in_command_file_number(in_command_file, "length");
	if (lg_name_phylip==NOT_IN_FILE)
		{
		printf("length of sequence names (number of characters) [10]. Type 0 if names have different lengths (separated from sequences by a space): ");
		lg_name_phylip=Answer(10);

		if (out_command_file) fprintf(out_command_file, "#specific to PHYLIP files\nname length: %d\n", lg_name_phylip);
		}
	}

n_lineages=Read_in_command_file_number(in_command_file, "lineages");
while (n_lineages<2)	/* NOT_IN_FILE is negative */
	{
	printf("number of lineages to compare [2]: ");
	n_lineages=Answer(2);
	}
if (out_command_file) fprintf(out_command_file, "lineages: %d\n", n_lineages);

lineage=Allocation(n_lineages+1, sizeof(struct sequences*), "lineage");
nb=Allocation(n_lineages+1, sizeof(int), "nb");
lineage_name=Allocation(n_lineages+1, sizeof(char*), "lineage_name");

for (i=0; i<n_lineages+1; i++) lineage[i]=Allocation(1, sizeof(struct sequences), "lineage[i]");

printf("For each sequence, give its lineage number (1 to %d), or 0 for the Outgroup lineage, or -1 to exclude a sequence from the computation:\n", n_lineages);

i=0;

if (out_command_file) fprintf(out_command_file, "\n#lineage affectation of sequences\n");

switch (file_format)
	{
	case mase:
		n=Nbseq_mase(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_mase(in, n, read, 0);
		break;
	case fasta:
		n=Nbseq_fasta_gde_phyltest(in, '>');
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_fasta_gde(in, n, read, '>');
		break;
	case gde:
		n=Nbseq_fasta_gde_phyltest(in, '%');
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_fasta_gde(in, n, read, '%');
		break;
	case phylip:
		n=Nbseq_phylip(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_phylip(in, n, 0, read, lg_name_phylip);
		break;
	case clustal:
		n=Nbseq_clustal(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_clustal(in, n, 0, read);
		break;
	case nexus:
		n=Nbseq_nexus(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		Read_nexus(in, n, 0, read);
		break;
	case lintre:
		n=Nbseq_lintre(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		the_end=Read_lintre(in, n, read);
		break;
	case phyltest:
		n=Nbseq_fasta_gde_phyltest(in, '#');
		read=Allocation(n, sizeof(struct sequences), "read");
		the_end=Read_phyltest(in, n, 0, read);
		break;
	case mega:
		n=Nbseq_mega(in);
		read=Allocation(n, sizeof(struct sequences), "read");
		the_end=Read_mega(in, n, 0, read);
		break;
	default:
		Leave_program("sorry, file format not implemented.");
	}

for (i=0; i<n; i++)
	{
	for (j=0; j<read[i].lg; j++) read[i].seq[j]=toupper(read[i].seq[j]);

	status=Read_in_command_file_number(in_command_file, read[i].name);

	while ((status<-1) || (status>n_lineages))	/* NOT_IN_FILE is <-1 */
		{
		printf("%s: ", read[i].name);
		Read_line(line, MAXLEN, stdin);
		if (line[0]) sscanf(line, "%d", &status);
		else status=-1;
		}

	if (out_command_file) fprintf(out_command_file, "%s: %d\n", read[i].name, status);

	if (status!=-1)
		{
		lineage[status][nb[status]]=Copy_sequence(read[i]);
		Free_sequences(&read[i], 1);
		nb[status]++;
		
		p=lineage[status];
		lineage[status]=Allocation(nb[status]+1, sizeof(struct sequences), "lineage[status]");

		for (j=0; j<nb[status]; j++)
			{
			lineage[status][j].lg=p[j].lg;
			lineage[status][j].name=p[j].name;
			lineage[status][j].com=p[j].com;
			lineage[status][j].seq=p[j].seq;
			}

		p=De_allocation(p);
		}
	}

read=De_allocation(read);

choices=Allocation(1, sizeof(struct options), "choices");
choices->first=1;
choices->several_sequences=0;

for (i=0; i<n_lineages+1; i++)
	if (!nb[i])
		{
		temp_string=Allocation(40, sizeof(char), "temp_string");
		sprintf(temp_string, "There are no sequences in lineage %d", i);
		Leave_program(temp_string);
		}
	else
		if (nb[i]>1) choices->several_sequences=1;

temp_string=Allocation(10, sizeof(char), "temp_string");

if (out_command_file) fprintf(out_command_file, "\n#lineage names\n");

for (i=0; i<n_lineages+1; i++)
	{
	if (i)
		{
		sprintf(temp_string, "lineage%d\0", i);
		lineage_name[i]=Read_in_command_file_string(in_command_file, temp_string);

		if (lineage_name[i]==NULL)
			{
			if (nb[i]>1) printf("name of lineage #%d [%s]: ", i, temp_string);
			else printf("name of lineage #%d [%s]: ", i, lineage[i][0].name);
			Read_line(line, MAXLEN, stdin);

			if (line[0]=='\n' || line[0]==0)
				{
				if (nb[i]>1)
					{
					lineage_name[i]=Allocation(10, sizeof(char), "lineage_name[i]");
					strcpy(lineage_name[i], temp_string);
					}
				else
					{
					lineage_name[i]=Allocation(strlen(lineage[i][0].name)+1, sizeof(char), "lineage_name[i]");
					strcpy(lineage_name[i], lineage[i][0].name);
					}
				}
			else
				{
				lineage_name[i]=Allocation(strlen(line), sizeof(char), "lineage_name[i]");
				strcpy(lineage_name[i], line);
				}
			}

		if (out_command_file) fprintf(out_command_file, "lineage%d: %s\n", i, lineage_name[i]);
		}
	else
		{
		lineage_name[i]=Read_in_command_file_string(in_command_file, "outgroup");
		if (lineage_name[i]==NULL)
			{
			printf("name of outgroup lineage [Outgroup]: ");
			Read_line(line, MAXLEN, stdin);
			if (line[0]=='\n' || line[0]==0)
				lineage_name[i]="Outgroup\0";
			else
				{
				lineage_name[i]=Allocation(strlen(line), sizeof(char), "lineage_name[i]");
				strcpy(lineage_name[i], line);
				}
			}

		if (out_command_file) fprintf(out_command_file, "outgroup: %s\n", lineage_name[i]);
		}
	}

temp_string=De_allocation(temp_string);

What_sequence_type(lineage[0][0].seq, choices, in_command_file);

printf("sequence type is assumed to be ");
switch (choices->seqtype)
	{
	case CDS:
		printf("coding DNA sequences (CDS)\n");
		if (out_command_file)
			{
			fprintf(out_command_file, "\nsequence type: cds\n");
			if (choices->code_mt) fprintf(out_command_file, "code: mitochondrial\n");
			else fprintf(out_command_file, "code: nuclear\n");
			}
		break;
	case NC:
		printf("non coding DNA sequences (NC)\n");
		if (out_command_file) fprintf(out_command_file, "sequence type: nc\n");
		break;
	case PROT:
		printf("protein sequences (PROT)\n");
		if (out_command_file) fprintf(out_command_file, "sequence type: prot\n");
	}

if (choices->seqtype!=PROT)
	{
	for (i=0; i<n_lineages+1; i++)
		for (j=0; j<nb[i]; j++)
			for (k=0; k<lineage[i][j].lg; k++)
				if (lineage[i][j].seq[k]=='U') lineage[i][j].seq[k]='T';
	}

if (choices->several_sequences)
	{
	choices->topology=Read_in_command_file_number(in_command_file, "topology");

	if (choices->topology==NOT_IN_FILE)
		{
		printf("take topology into account (you need a tree file) ([1]/0) ? ");
		choices->topology=Answer(1);
		}
	}
else
	choices->topology=0;

if (out_command_file) fprintf(out_command_file, "topology: %d\n", choices->topology);

if (choices->topology)
	{
	read_name=Read_in_command_file_string(in_command_file, "tree");

	if (read_name==NULL)
		{
		i=0;
		while (name_file[i] && name_file[i]!='.')
			{
			line[i]=name_file[i];
			i++;
			}
		line[i]=0;

		read_name=Allocation(i+5, sizeof(char), "read_name");
		if (WINDOWS)	sprintf(read_name, "%s.tre", line);
		else	sprintf(read_name, "%s.tree", line);

		printf("NEWICK file containing a rooted tree [%s]: ", read_name);
		Read_line(line, MAXLEN, stdin);
		if (line[0] && line[0]!='\n')	sscanf(line, "%s", read_name);
		}

	if (out_command_file) fprintf(out_command_file, "tree file: %s\n", read_name);

	ftree=Open_file("r", read_name);

	choices->support_limit=Read_in_command_file_number(in_command_file, "threshold");

	if (choices->support_limit==NOT_IN_FILE)
		{
		printf("ignore lowly supported branches (1/[0]) ? ");
		choices->support_limit=Answer(0);

		if (choices->support_limit)
			{
			printf("threshold [50]:");
			choices->support_limit=Answer(50);
			}
		}
	
	tree_parenthesis=Read_tree(ftree, choices);
	}
else
	{
	ftree=NULL;
	choices->support_limit=0;
	}

if (out_command_file) fprintf(out_command_file, "threshold: %d\n", choices->support_limit);

read_name=Read_in_command_file_string(in_command_file, "text");
if (read_name)
	{
	if (read_name[0]=='0' && read_name[1]==0) out_text=stdout;
	else out_text=Open_file("a", read_name);
	}
else
	{
	read_name=Allocation(NAMELEN, sizeof(char), "read_name");
	printf("output text file [screen]: ");
	gets(read_name);
	if (read_name[0]=='\n' || read_name[0]==0) out_text=stdout;
	else out_text=Open_file("a", read_name);
	}

if (out_command_file)
	{
	if (out_text==stdout) fprintf(out_command_file, "out text file: 0\n");
	else fprintf(out_command_file, "out text file: %s\n", read_name);
	}

read_name=De_allocation(read_name);

read_name=Read_in_command_file_string(in_command_file, "table");
if (read_name)
	{
	if (read_name[0]=='0' && read_name[1]==0) out_table=NULL;	/* the command file can specify: no table output by giving it the name "0" */
	else out_table=Open_file("w", read_name);
	}
else
	{
	read_name=Allocation(NAMELEN, sizeof(char), "read_name");
	printf("output table file [no such file]: ");
	gets(read_name);
	if (read_name[0]=='\n' || read_name[0]==0) out_table=NULL;
	else out_table=Open_file("w", read_name);
	}

if (out_command_file)
	{
	if (out_table==NULL) fprintf(out_command_file, "out table file: 0\n");
	else fprintf(out_command_file, "out table file: %s\n", read_name);
	}

read_name=De_allocation(read_name);

switch (choices->seqtype)
	{
	case CDS:
		choices->compute_cds[ks]=Read_in_command_file_number(in_command_file, "ks");
		if (choices->compute_cds[ks]==NOT_IN_FILE)
			{
			printf("compare synonymous substitutions = Ks ([1]/0) ? ");
			choices->compute_cds[ks]=Answer(1);
			}
		if (out_command_file) fprintf(out_command_file, "ks: %d\n", choices->compute_cds[ks]);

		choices->compute_cds[ka]=Read_in_command_file_number(in_command_file, "ka");
		if (choices->compute_cds[ka]==NOT_IN_FILE)
			{
			printf("compare non synonymous substitutions = Ka ([1]/0) ? ");
			choices->compute_cds[ka]=Answer(1);
			}
		if (out_command_file) fprintf(out_command_file, "ka: %d\n", choices->compute_cds[ka]);

		choices->compute_cds[as]=Read_in_command_file_number(in_command_file, "as");
		if (choices->compute_cds[as]==NOT_IN_FILE)
			{
			printf("compare synonymous substitutions = As (1/[0]) ? ");
			choices->compute_cds[as]=Answer(0);
			}
		if (out_command_file) fprintf(out_command_file, "as: %d\n", choices->compute_cds[as]);

		choices->compute_cds[ba]=Read_in_command_file_number(in_command_file, "ba");
		if (choices->compute_cds[ba]==NOT_IN_FILE)
			{
			printf("compare synonymous substitutions = Ba (1/[0]) ? ");
			choices->compute_cds[ba]=Answer(0);
			}
		if (out_command_file) fprintf(out_command_file, "ba: %d\n", choices->compute_cds[ba]);

		choices->compute_cds[b4]=Read_in_command_file_number(in_command_file, "b4");
		if (choices->compute_cds[b4]==NOT_IN_FILE)
			{
			printf("compare synonymous substitutions = B4 (1/[0]) ? ");
			choices->compute_cds[b4]=Answer(0);
			}
		if (out_command_file) fprintf(out_command_file, "b4: %d\n", choices->compute_cds[b4]);

		break;
	case NC:
		read_name=Read_in_command_file_string(in_command_file, "distance");
		if (read_name)
			{
			if (read_name[0]=='j' && read_name[1]=='c') choices->compute_nc=jc;
			else choices->compute_nc=k2;
			}
		else
			{
			printf("distance method:\n");
			printf("Jukes-Cantor 1 parameter (%d), or Kimura 2-parameter [%d]? ", jc, k2);
			choices->compute_nc=Answer(k2);
			}

		if (out_command_file)
			{
			fprintf(out_command_file, "distance: ");
			switch (choices->compute_nc)
				{
				case jc:
					fprintf(out_command_file, "jc\n");
					break;
				case k2:
					fprintf(out_command_file, "k2\n");
					break;
				default:
					fprintf(out_command_file, "\n");
				}
			}
	}

choices->verbose=Read_in_command_file_number(in_command_file, "verbose");
if (choices->verbose==NOT_IN_FILE)
	{
	printf("Do you want to print all computed variances ([0]/1) ? ");
	choices->verbose=Answer(0);
	if (out_command_file) fprintf(out_command_file, "verbose: %d\n", choices->verbose);
	}

switch (choices->seqtype)
	{
	case CDS:
	case NC:
		Allgap(lineage, nb, n_lineages, bases);
		break;
	case PROT:
		Allgap(lineage, nb, n_lineages, aminoacids);
	}

for (i=0; i<n_lineages+1; i++)
	for (j=0; j<nb[i]; j++)
		switch (choices->seqtype)
			{
			case CDS:
				lineage[i][j]=Nogap_cds(lineage[i][j], bases);
				break;
			case NC:
				lineage[i][j]=Nogap(lineage[i][j], bases);
				break;
			case PROT:
				lineage[i][j]=Nogap(lineage[i][j], aminoacids);
			}

num_lin[0]=0;

fprintf(out_text, "\n\n## treatement of %d lineages from file %s ##\n\n", n_lineages, name_file);

if (out_table)
	{
	fprintf(out_table, "outgroup\tnb_seq_outgroup\tlineage1\tnb_seq_lin1\tlineage2\tnb_seq_lin2\tweighting\tsupport_limit\tseq_type\t");

	if (choices->seqtype<PROT)
		for (i=0; i<3; i++)
			{
			fprintf(out_table, "GC_tot_%d\t", i);
			if (choices->seqtype==CDS) for (j=1; j<=3; j++) fprintf(out_table, "GC%d_%d\t", j, i);
			}

	if (choices->seqtype==CDS)
		{
		if (choices->compute_cds[ks]) fprintf(out_table, "nb_sites_Ks\tKs1\tKs2\tdKs\tsd_dKs\tratio_Ks\tP_Ks\t");
		if (choices->compute_cds[as]) fprintf(out_table, "nb_sites_As\tAs1\tAs2\tdAs\tsd_dAs\tratio_As\tP_As\t");
		if (choices->compute_cds[b4]) fprintf(out_table, "nb_sites_B4\tB41\tB42\tdB4\tsd_dB4\tratio_B4\tP_B4\t");
		if (choices->compute_cds[ka]) fprintf(out_table, "nb_sites_Ka\tKa1\tKa2\tdKa\tsd_dKa\tratio_Ka\tP_Ka\t");
		if (choices->compute_cds[ba]) fprintf(out_table, "nb_sites_Ba\tBa1\tBa2\tdBa\tsd_dBa\tratio_Ba\tP_Ba\t");
		fprintf(out_table, "\n");
		}
	else
		fprintf(out_table, "nb_sites\tK1\tK2\tdK\tsd_dK\tratio\tP\t\n");
	}

for (num_lin[1]=1; num_lin[1]<n_lineages; num_lin[1]++)
	for (num_lin[2]=num_lin[1]+1; num_lin[2]<n_lineages+1; num_lin[2]++)
		if (choices->seqtype==CDS) Treatement_cds(out_text, out_table, tree_parenthesis, lineage, nb, num_lin, choices, lineage_name);
		else Treatement_nc(out_text, out_table, tree_parenthesis, lineage, nb, num_lin, choices, lineage_name);

printf("\nNormal end of the program.\n");

if (WINDOWS)
	{
	printf("type enter\n");
	Read_line(line, MAXLEN, stdin);
	}

fclose(in);
fclose(out_text);
if (out_table) fclose(out_table);
}

/*****************************************************************************/

void Treatement_cds(FILE *out_text, FILE *out_table, char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin, struct options *choices, char ** lineage_name)
{
struct sequences to_treat[2];
int i, j, k, r, rr;
double l;
struct parameters_lwl ***para;	/*results of Ks-Ka computations for all pairs lineage - outgroup*/
struct parameters_lwl p12, p34;	/*computation intermediaries for covariance*/
double cov;	/*computation intermediaries for covariance*/
double Sls[2], Sla[2], Sl4[2]; /*length sums*/
double Sllcov_ks, Sllcov_ka;	/*computation intermediaries for variance of the mean difference*/
double Kst[2], Kat[2];
double dKs=1./PRECISION, var_dKs=0, dKa=1./PRECISION, var_dKa=0;	/*mean differences of Ks-Ka and their variances*/
double **weight;
double weight_pair;
int pb_ks=1-choices->compute_cds[ks], pb_ka=1-choices->compute_cds[ka], pb_as=1-choices->compute_cds[as], pb_ba=1-choices->compute_cds[ba], pb_b4=1-choices->compute_cds[b4];
double Sllcov_as, Ast[2], dAs, var_dAs=0;	/*comparison of synonymous transitions*/
double Sllcov_ba, Bat[2], dBa, var_dBa=0;	/*comparison of non synonymous transversions*/
double Sllcov_b4, B4t[2], dB4, var_dB4=0;	/*comparison of synonymous transversions*/
double gc_mean[3][4], *gc;
char one_char;

if (choices->several_sequences)
	{
	if (choices->topology)
		{
		weight=Weighting(tree_parenthesis, lineage, nb, num_lin);

		fprintf(out_text, "\ntopologic weight of sequences:\n");
		}
	else
		{
		weight=Allocation(3, sizeof(double*), "weight");
		for (i=0; i<3; i++)
			{
			weight[i]=Allocation(nb[num_lin[i]], sizeof(double), "weight[i]");
			for (j=0; j<nb[num_lin[i]]; j++) weight[i][j]=Divide(1., (double)nb[num_lin[i]]);
			}

		fprintf(out_text, "equally weighted sequences:\n");
		}

	for (i=0; i<3; i++)
		{
		fprintf(out_text, "> %s:\n", lineage_name[num_lin[i]]);

		for (j=0; j<nb[num_lin[i]]; j++)
			fprintf(out_text, "%s\t%g\n", lineage[num_lin[i]][j].name, weight[i][j]);

		if (out_table) fprintf(out_table, "%s\t%d\t", lineage_name[num_lin[i]], nb[num_lin[i]]);
		}

	fprintf(out_text, "\n");
	}
else
	{
	weight=Allocation(3, sizeof(double*), "weight");
	for (i=0; i<3; i++)
		{
		weight[i]=Allocation(nb[num_lin[i]], sizeof(double), "weight[i]");
		for (j=0; j<nb[num_lin[i]]; j++) weight[i][j]=Divide(1., (double)nb[num_lin[i]]);
		}

	if (out_table) for (i=0; i<3; i++) fprintf(out_table, "%s\t%d\t", lineage_name[num_lin[i]], nb[num_lin[i]]);
	}

if (out_table)
	{
	fprintf(out_table, "%d\t%d\t", choices->topology, choices->support_limit);
	switch (choices->seqtype)
		{
		case NC:
			fprintf(out_table, "NC\t");
			break;
		case CDS:
			fprintf(out_table, "CDS\t");
			break;
		case PROT:
			fprintf(out_table, "PROT\t");
		}
	}

for (i=0; i<3; i++)
	for (j=0; j<4; j++) gc_mean[i][j]=0;

for (r=0; r<nb[0]; r++)
	{
	gc=GC(lineage[0][r], choices);
	for (i=0; i<4; i++) gc_mean[0][i]+=weight[0][r]*gc[i];
	gc=De_allocation(gc);
	}

/*all comparisons outgroup sequences - other sequences*/

para=Allocation(2, sizeof(struct parameters_lwl**), "para");

for (i=0; i<2; i++)	/*two lineages*/
	{

/*first time through to count things and compute length sums*/

	para[i]=Allocation(nb[0], sizeof(struct parameters_lwl*), "para[i]");
	
	Sls[i]=Sla[i]=Sl4[i]=0;
	
	for (r=0; r<nb[0]; r++)		/*nb[0] outgroup sequences*/
		{
		para[i][r]=Allocation(nb[num_lin[i+1]], sizeof(struct parameters_lwl), "para[i][r]");

		to_treat[0]=lineage[0][r];

		for (j=0; j<nb[num_lin[i+1]]; j++)	/*nb[] sequences per lineage*/
			{
			to_treat[1]=lineage[num_lin[i+1]][j];
			para[i][r][j]=Lwl(to_treat, choices);

			if (choices->verbose)
				{
				if (choices->compute_cds[ks])
					{
					if (para[i][r][j].ks>=0) fprintf(out_text, "Ks (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].ks, sqrt(para[i][r][j].vks));
					else fprintf(out_text, "Ks (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
					}
				if (choices->compute_cds[ka])
					{
					if (para[i][r][j].ka>=0) fprintf(out_text, "Ka (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].ka, sqrt(para[i][r][j].vka));
					else fprintf(out_text, "Ka (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
					}
				if (choices->compute_cds[as])
					{
					if (para[i][r][j].as>=0) fprintf(out_text, "As (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].as, sqrt(para[i][r][j].vas));
					else fprintf(out_text, "As (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
					}
				if (choices->compute_cds[ba])
					{
					if (para[i][r][j].ba>=0) fprintf(out_text, "Ba (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].ba, sqrt(para[i][r][j].vba));
					else fprintf(out_text, "Ba (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
					}
				if (choices->compute_cds[b4])
					{
					if (para[i][r][j].b4>=0) fprintf(out_text, "B4 (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].b4, sqrt(para[i][r][j].vb4));
					else fprintf(out_text, "B4 (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
					}
				}

			if (choices->first) choices->first=0;
			if (para[i][r][j].ks<0) pb_ks=1;
			if (para[i][r][j].ka<0) pb_ka=1;
			if (para[i][r][j].as<0) pb_as=1;
			if (para[i][r][j].ba<0) pb_ba=1;
			if (para[i][r][j].b4<0) pb_b4=1;

			weight_pair=weight[0][r]*weight[i+1][j];

			if (!pb_ks || !pb_as) Sls[i]+=weight_pair*para[i][r][j].ls;
			if (!pb_ka || !pb_ba) Sla[i]+=weight_pair*para[i][r][j].la;
			if (!pb_b4) Sl4[i]+=weight_pair*para[i][r][j].l4;

			gc=GC(to_treat[1], choices);
			for (k=0; k<4; k++) gc_mean[i+1][k]+=weight_pair*gc[k];
			gc=De_allocation(gc);
			}
		}

/*second time through to compute weighted means and variances*/

	Kst[i]=Kat[i]=0;
	Ast[i]=Bat[i]=B4t[i]=0.;

	for (r=0; r<nb[0]; r++)	/*nb[0] outgroup sequences*/
		for (j=0; j<nb[num_lin[i+1]]; j++)	/*nb[] sequences per lineage*/
			{
			weight_pair=Divide(para[i][r][j].ls, Sls[i]);
			weight_pair*=weight[0][r]*weight[i+1][j];

			if (!pb_ks) Kst[i]+=weight_pair*para[i][r][j].ks;
			if (!pb_as) Ast[i]+=weight_pair*para[i][r][j].as;

			weight_pair=Divide(para[i][r][j].la, Sla[i]);
			weight_pair*=weight[0][r]*weight[i+1][j];

			if (!pb_ka) Kat[i]+=weight_pair*para[i][r][j].ka;
			if (!pb_ba) Bat[i]+=weight_pair*para[i][r][j].ba;

			weight_pair=Divide(para[i][r][j].l4, Sl4[i]);
			weight_pair*=weight[0][r]*weight[i+1][j];

			if (!pb_b4) B4t[i]+=weight_pair*para[i][r][j].b4;
			}

	for (j=0; j<nb[num_lin[i+1]]; j++)
		for (k=0; k<nb[num_lin[i+1]]; k++)
			{
			to_treat[0]=lineage[num_lin[i+1]][j];
			to_treat[1]=lineage[num_lin[i+1]][k];
			p12=Lwl(to_treat, choices);
			if (p12.ks<0) pb_ks=1;
			if (p12.ka<0) pb_ka=1;
			if (p12.as<0) pb_as=1;
			if (p12.ba<0) pb_ba=1;
			if (p12.b4<0) pb_b4=1;

			for (r=0; r<nb[0]; r++)
				for (rr=0; rr<nb[0]; rr++)
						{
						to_treat[0]=lineage[0][r];
						to_treat[1]=lineage[0][rr];
						p34=Lwl(to_treat, choices);
						if (p34.ks<0) pb_ks=1;
						if (p34.ka<0) pb_ka=1;
						if (p34.as<0) pb_as=1;
						if (p34.ba<0) pb_ba=1;
						if (p34.b4<0) pb_b4=1;
						weight_pair=Divide(para[i][r][j].ls, Sls[i]);
						weight_pair*=Divide(para[i][rr][k].ls, Sls[i]);
						weight_pair*=weight[0][r]*weight[0][rr]*weight[i+1][j]*weight[i+1][k];

						if (!pb_ks)
							{
							cov=CovKs(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34);
							if (choices->verbose && (j!=k || r!=rr))
								{
								if (cov>=0) fprintf(out_text, "Ks covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
								else fprintf(out_text, "Ks covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
								}
							var_dKs+=weight_pair*cov;
							}

						if (!pb_as)
							{
							cov=CovAs(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34);
							if (choices->verbose && (j!=k || r!=rr))
								{
								if (cov>=0) fprintf(out_text, "As covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
								else fprintf(out_text, "As covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
								}
							var_dAs+=weight_pair*cov;
							}

						weight_pair=Divide(para[i][r][j].la, Sla[i]);
						weight_pair*=Divide(para[i][rr][k].la, Sla[i]);
						weight_pair*=weight[0][r]*weight[0][rr]*weight[i+1][j]*weight[i+1][k];

						if (!pb_ka)
							{
							cov=CovKa(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34);
							if (choices->verbose && (j!=k || r!=rr))
								{
								if (cov>=0) fprintf(out_text, "Ka covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
								else fprintf(out_text, "Ka covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
								}
							var_dKa+=weight_pair*cov;
							}

						if (!pb_ba)
							{
							cov=CovBa(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34);
							if (choices->verbose && (j!=k || r!=rr))
								{
								if (cov>=0) fprintf(out_text, "Ba covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
								else fprintf(out_text, "Ba covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
								}
							var_dBa+=weight_pair*cov;
							}

						weight_pair=Divide(para[i][r][j].l4, Sl4[i]);
						weight_pair*=Divide(para[i][rr][k].l4, Sl4[i]);
						weight_pair*=weight[0][r]*weight[0][rr]*weight[i+1][j]*weight[i+1][k];

						if (!pb_b4)
							{
							cov=CovB4(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34);
							if (choices->verbose && (j!=k || r!=rr))
								{
								if (cov>=0) fprintf(out_text, "B4 covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
								else fprintf(out_text, "B4 covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
								}
							var_dB4+=weight_pair*cov;
							}
						}
			}
	}

Sllcov_ks=0;
Sllcov_ka=0;
Sllcov_as=0;
Sllcov_ba=0;
Sllcov_b4=0;

for (i=0; i<nb[num_lin[1]]; i++)
	for (j=0; j<nb[num_lin[2]]; j++)
		{
		to_treat[0]=lineage[num_lin[1]][i];
		to_treat[1]=lineage[num_lin[2]][j];
		p12=Lwl(to_treat, choices);
		if (p12.ks<0) pb_ks=1;
		if (p12.ka<0) pb_ka=1;
		if (p12.as<0) pb_as=1;
		if (p12.ba<0) pb_ba=1;
		if (p12.b4<0) pb_b4=1;

		for (r=0; r<nb[0]; r++)
			for (rr=0; rr<nb[0]; rr++)
				{
				to_treat[0]=lineage[0][r];
				to_treat[1]=lineage[0][rr];
				p34=Lwl(to_treat, choices);
				if (p34.ks<0) pb_ks=1;
				if (p34.ka<0) pb_ka=1;
				if (p34.as<0) pb_as=1;
				if (p34.ba<0) pb_ba=1;
				if (p34.b4<0) pb_b4=1;

				weight_pair=Divide(para[0][r][i].ls, Sls[0]);
				weight_pair*=Divide(para[1][rr][j].ls, Sls[1]);
				weight_pair*=weight[0][r]*weight[0][rr]*weight[1][i]*weight[2][j];

				if (!pb_ks)
					{
					cov=CovKs(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34);
					Sllcov_ks+=weight_pair*cov;
					}

				if (!pb_as)
					{
					cov=CovAs(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34);
					Sllcov_as+=weight_pair*cov;
					}

				weight_pair=Divide(para[0][r][i].la, Sla[0]);
				weight_pair*=Divide(para[1][rr][j].la, Sla[1]);
				weight_pair*=weight[0][r]*weight[0][rr]*weight[1][i]*weight[2][j];

				if (!pb_ka)
					{
					cov=CovKa(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34);
					Sllcov_ka+=weight_pair*cov;
					}

				if (!pb_ba)
					{
					cov=CovBa(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34);
					Sllcov_ba+=weight_pair*cov;
					}

				weight_pair=Divide(para[0][r][i].l4, Sl4[0]);
				weight_pair*=Divide(para[1][rr][j].l4, Sl4[1]);
				weight_pair*=weight[0][r]*weight[0][rr]*weight[1][i]*weight[2][j];

				if (!pb_b4)
					{
					cov=CovB4(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34);
					Sllcov_b4+=weight_pair*cov;
					}
				}
		}

for (i=0; i<3; i++) weight[i]=De_allocation(weight[i]);
weight=De_allocation(weight);

if (!pb_ks)
	{
	if (choices->verbose) fprintf(out_text, "Ks covariance between ingroups=%g (sd=%g)\n", Sllcov_ks, sqrt(Sllcov_ks));
	dKs=Kst[0] - Kst[1];
	if (dKs==0) pb_ks=1;
	}
if (!pb_ka)
	{
	if (choices->verbose) fprintf(out_text, "Ka covariance between ingroups=%g (sd=%g)\n", Sllcov_ka, sqrt(Sllcov_ka));
	dKa=Kat[0] - Kat[1];
	if (dKa==0) pb_ka=1;
	}
if (!pb_as)
	{
	if (choices->verbose) fprintf(out_text, "As covariance between ingroups=%g (sd=%g)\n", Sllcov_as, sqrt(Sllcov_as));
	dAs=Ast[0] - Ast[1];
	if (dAs==0) pb_as=1;
	}
if (!pb_ba)
	{
	if (choices->verbose) fprintf(out_text, "Ba covariance between ingroups=%g (sd=%g)\n", Sllcov_ba, sqrt(Sllcov_ba));
	dBa=Bat[0] - Bat[1];
	if (dBa==0) pb_ba=1;
	}
if (!pb_b4)
	{
	if (choices->verbose) fprintf(out_text, "B4 covariance between ingroups=%g (sd=%g)\n", Sllcov_b4, sqrt(Sllcov_b4));
	dB4=B4t[0] - B4t[1];
	if (dB4==0) pb_b4=1;
	}

for (i=0; i<3; i++)
	{
	fprintf(out_text, "mean GC ");
	fprintf(out_text, "%s: ", lineage_name[num_lin[i]]);
	fprintf(out_text, "all positions=%.3lf ", gc_mean[i][0]);
	for (j=1; j<4; j++)
		fprintf(out_text, "position%d=%.3lf ", j, gc_mean[i][j]);
	fprintf(out_text, "\n");

	if (out_table)
		for (j=0; j<4; j++)
			fprintf(out_table, "%g\t", gc_mean[i][j]);
	}
fprintf(out_text, "\n");

if (!pb_ks)
	{
	var_dKs-=2*Sllcov_ks;
	if (var_dKs<=0) pb_ks=1;
	}

if (!pb_ks)
	{
	fprintf(out_text, "mean number of synonymous sites compared: %.1lf\n", Divide(Sls[0]+Sls[1], 2.));
	for (i=0; i<2; i++) fprintf(out_text, "mean Ks %s: %g\n", lineage_name[num_lin[i+1]], Kst[i]);

	fprintf(out_text, "dKs: %g\tstandard deviation (sd): %g\tratio dKs/sd: %g\n", dKs, sqrt(var_dKs), Divide(dKs, sqrt(var_dKs)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dKs, sqrt(var_dKs))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", Divide(Sls[0]+Sls[1], 2.), Kst[0], Kst[1], dKs, sqrt(var_dKs), Divide(dKs, sqrt(var_dKs)), Exact_probability(Divide(dKs, sqrt(var_dKs))));
	}
else
	{
	if (choices->compute_cds[ks])
		{
		if (dKs==0)
			{
			fprintf(out_text, "mean number of synonymous sites compared: %.1lf\n", Divide(Sls[0]+Sls[1], 2.));
			for (i=0; i<2; i++) fprintf(out_text, "mean Ks %s: %g\n", lineage_name[num_lin[i+1]], Kst[i]);
			fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

			if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", Divide(Sls[0]+Sls[1], 2.), Kst[0], Kst[1]);
			}
		else
			{
			fprintf(out_text, "** impossible to compute Ks. saturation? **\n\n");
			if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
			}
		}
	}

if (!pb_as)
	{
	var_dAs-=2*Sllcov_as;
	if (var_dAs<=0) pb_as=1;
	}

if (!pb_as)
	{
	fprintf(out_text, "mean number of synonymous sites compared: %.1lf\n", Divide(Sls[0]+Sls[1], 2.));
	for (i=0; i<2; i++) fprintf(out_text, "mean As %s: %g\n", lineage_name[num_lin[i+1]], Ast[i]);

	fprintf(out_text, "dAs: %g\tstandard deviation (sd): %g\tratio dAs/sd: %g\n", dAs, sqrt(var_dAs), Divide(dAs, sqrt(var_dAs)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dAs, sqrt(var_dAs))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", Divide(Sls[0]+Sls[1], 2.), Ast[0], Ast[1], dAs, sqrt(var_dAs), Divide(dAs, sqrt(var_dAs)), Exact_probability(Divide(dAs, sqrt(var_dAs))));
	}
else
	{
	if (choices->compute_cds[as])
		{
		if (dAs==0)
			{
			fprintf(out_text, "mean number of synonymous sites compared: %.1lf\n", Divide((Sls[0]+Sls[1]), 2.));
			for (i=0; i<2; i++) fprintf(out_text, "mean As %s: %g\n", lineage_name[num_lin[i+1]], Ast[i]);
			fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

			if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", Divide(Sls[0]+Sls[1], 2.), Ast[0], Ast[1]);
			}
		else
			{
			fprintf(out_text, "** impossible to compute As. saturation? **\n\n");
			if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
			}
		}
	}

if (!pb_b4)
	{
	var_dB4-=2*Sllcov_b4;
	if (var_dB4<=0) pb_b4=1;
	}

if (!pb_b4)
	{
	fprintf(out_text, "mean number of totally degenerate sites compared: %.1lf\n", Divide(Sl4[0]+Sl4[1], 2.));
	for (i=0; i<2; i++) fprintf(out_text, "mean B4 %s: %g\n", lineage_name[num_lin[i+1]], B4t[i]);

	fprintf(out_text, "dB4: %g\tstandard deviation (sd): %g\tratio dB4/sd: %g\n", dB4, sqrt(var_dB4), Divide(dB4, sqrt(var_dB4)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dB4, sqrt(var_dB4))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", Divide((Sl4[0]+Sl4[1]), 2.), B4t[0], B4t[1], dB4, sqrt(var_dB4), Divide(dB4, sqrt(var_dB4)), Exact_probability(Divide(dB4, sqrt(var_dB4))));
	}
else
	{
	if (choices->compute_cds[b4])
		{
		if (dB4==0)
			{
			fprintf(out_text, "mean number of totally degenerate sites compared: %.1lf\n", Divide((Sl4[0]+Sl4[1]), 2.));
			for (i=0; i<2; i++) fprintf(out_text, "mean B4 %s: %g\n", lineage_name[num_lin[i+1]], B4t[i]);
			fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

			if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", Divide(Sl4[0]+Sl4[1], 2.), B4t[0], B4t[1]);
			}
		else
			{
			fprintf(out_text, "** impossible to compute B4. saturation? **\n\n");
			if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
			}
		}
	}

if (!pb_ka)
	{
	var_dKa-=2*Sllcov_ka;
	if (var_dKa<=0) pb_ka=1;
	}

if (!pb_ka)
	{
	fprintf(out_text, "mean number of non synonymous sites compared: %.1lf\n", Divide(Sla[0]+Sla[1], 2.));
	for (i=0; i<2; i++) fprintf(out_text, "mean Ka %s: %g\n", lineage_name[num_lin[i+1]], Kat[i]);

	fprintf(out_text, "dKa: %g\tstandard deviation (sd): %g\tratio dKa/sd: %g\n", dKa, sqrt(var_dKa), Divide(dKa, sqrt(var_dKa)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dKa, sqrt(var_dKa))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", Divide(Sla[0]+Sla[1], 2.), Kat[0], Kat[1], dKa, sqrt(var_dKa), Divide(dKa, sqrt(var_dKa)), Exact_probability(Divide(dKa, sqrt(var_dKa))));
	}
else
	{
	if (choices->compute_cds[ka])
		{
		if (dKa==0)
			{
			fprintf(out_text, "mean number of non synonymous sites compared: %.1lf\n", Divide(Sla[0]+Sla[1], 2.));
			for (i=0; i<2; i++) fprintf(out_text, "mean Ka %s: %g\n", lineage_name[num_lin[i+1]], Kat[i]);
			fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

			if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", Divide(Sla[0]+Sla[1], 2.), Kat[0], Kat[1]);
			}
		else
			{
			fprintf(out_text, "** impossible to compute Ka. saturation? **\n\n");
			if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
			}
		}
	}

if (!pb_ba)
	{
	var_dBa-=2*Sllcov_ba;
	if (var_dBa<=0) pb_ba=1;
	}

if (!pb_ba)
	{
	fprintf(out_text, "mean number of non synonymous sites compared: %.1lf\n", Divide(Sla[0]+Sla[1], 2.));
	for (i=0; i<2; i++) fprintf(out_text, "mean Ba %s: %g\n", lineage_name[num_lin[i+1]], Bat[i]);

	fprintf(out_text, "dBa: %g\tstandard deviation (sd): %g\tratio dBa/sd: %g\n", dBa, sqrt(var_dBa), Divide(dBa, sqrt(var_dBa)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dBa, sqrt(var_dBa))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", Divide(Sla[0]+Sla[1], 2.), Bat[0], Bat[1], dBa, sqrt(var_dBa), Divide(dBa, sqrt(var_dBa)), Exact_probability(Divide(dBa, sqrt(var_dBa))));
	}
else
	{
	if (choices->compute_cds[ba])
		{
		if (dBa==0)
			{
			fprintf(out_text, "mean number of non synonymous sites compared: %.1lf\n", Divide(Sla[0]+Sla[1], 2.));
			for (i=0; i<2; i++) fprintf(out_text, "mean Ba %s: %g\n", lineage_name[num_lin[i+1]], Bat[i]);
			fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

			if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", Divide(Sla[0]+Sla[1], 2.), Bat[0], Bat[1]);
			}
		else
			{
			fprintf(out_text, "** impossible to compute Ba. saturation? **\n\n");
			if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
			}
		}
	}

for (i=0; i<2; i++)
	{
	for (j=0; j<nb[0]; j++) para[i][j]=De_allocation(para[i][j]);
	para[i]=De_allocation(para[i]);
	}
para=De_allocation(para);

if (out_table) fprintf(out_table, "\n");
}

/*****************************************************************************/

void Treatement_nc(FILE *out_text, FILE *out_table, char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin, struct options *choices, char ** lineage_name)
{
struct sequences to_treat[2];
int i, j, k, r, rr;
double l;
struct parameters_k ***para;	/*results of distance computation for all pairs lineages - outgroup*/
struct parameters_k p12, p34;	/*intermediairies for covariance computation*/
double cov;	/*intermediairy for covariance computation*/
double Sllcov;	/*intermediairy for computation of the variance of the mean difference*/
double Kt[2];
double dK, var_dK=0;	/*mean K difference and its variance*/
double **weight;
double weight_pair;
int pb_k=0;
double gc_mean[3], *gc;
char one_char;


if (choices->several_sequences)
	{
	if (choices->topology)
		{
		weight=Weighting(tree_parenthesis, lineage, nb, num_lin);

		fprintf(out_text, "\ntopologic weight of sequences:\n");
		}
	else
		{
		weight=Allocation(3, sizeof(double*), "weight");
		for (i=0; i<3; i++)
			{
			weight[i]=Allocation(nb[num_lin[i]], sizeof(double), "weight[i]");
			for (j=0; j<nb[num_lin[i]]; j++) weight[i][j]=Divide(1., (double)nb[num_lin[i]]);
			}

		fprintf(out_text, "equally weighted sequences:\n");
		}

	for (i=0; i<3; i++)
		{
		fprintf(out_text, "> %s:\n", lineage_name[num_lin[i]]);

		for (j=0; j<nb[num_lin[i]]; j++)
			fprintf(out_text, "%s\t%g\n", lineage[num_lin[i]][j].name, weight[i][j]);

		if (out_table) fprintf(out_table, "%s\t%d\t", lineage_name[num_lin[i]], nb[num_lin[i]]);
		}

	fprintf(out_text, "\n");
	}
else
	{
	weight=Allocation(3, sizeof(double*), "weight");
	for (i=0; i<3; i++)
		{
		weight[i]=Allocation(nb[num_lin[i]], sizeof(double), "weight[i]");
		for (j=0; j<nb[num_lin[i]]; j++) weight[i][j]=Divide(1., (double)nb[num_lin[i]]);
		}

	if (out_table) for (i=0; i<3; i++) fprintf(out_table, "%s\t%d\t", lineage_name[num_lin[i]], nb[num_lin[i]]);
	}

if (out_table)
	{
	fprintf(out_table, "%d\t%d\t", choices->topology, choices->support_limit);
	switch (choices->seqtype)
		{
		case NC:
			fprintf(out_table, "NC\t");
			break;
		case CDS:
			fprintf(out_table, "CDS\t");
			break;
		case PROT:
			fprintf(out_table, "PROT\t");
		}
	}

if (!choices->seqtype)
	{
	for (i=0; i<3; i++) gc_mean[i]=0;

	for (r=0; r<nb[0]; r++)
		{
		gc=GC(lineage[0][r], choices);
		gc_mean[0]+=weight[0][r]*gc[0];
		gc=De_allocation(gc);
		}
	}

/*all comparisons outgroup sequences - other sequences*/

para=Allocation(2, sizeof(struct parameters_k**), "para");

for (i=0; i<2; i++)	/*two lineages*/
	{

/*first time through to count things and compute length sums*/

	para[i]=Allocation(nb[0], sizeof(struct parameters_k*), "para[i]");

	Kt[i]=0;

	for (r=0; r<nb[0]; r++)		/*nb[0] outgroup sequences*/
		{
		para[i][r]=Allocation(nb[num_lin[i+1]], sizeof(struct parameters_k), "para[i][r]");

		to_treat[0]=lineage[0][r];

		for (j=0; j<nb[num_lin[i+1]]; j++)	/*nb[] sequences per lineage*/
			{
			to_treat[1]=lineage[num_lin[i+1]][j];
			if (choices->seqtype) para[i][r][j]=Poisson(to_treat, PROT_POISSON);
			else para[i][r][j]=NC_dist(to_treat, choices->compute_nc);
			if (para[i][r][j].k<0) pb_k=1;

			if (choices->verbose)
				{
				if (para[i][r][j].k>=0) fprintf(out_text, "K (%s, %s) = %g\tsd = %g\n", to_treat[0].name, to_treat[1].name, para[i][r][j].k, sqrt(para[i][r][j].vk));
				else fprintf(out_text, "K (%s, %s) not computable\n", to_treat[0].name, to_treat[1].name);
				}

			weight_pair=weight[0][r]*weight[i+1][j];

			if (!pb_k) Kt[i]+=weight_pair*para[i][r][j].k;

			if (!choices->seqtype)
				{
				gc=GC(to_treat[1], choices);
				gc_mean[i+1]+=weight_pair*gc[0];
				gc=De_allocation(gc);
				}
			}
		}

	for (j=0; j<nb[num_lin[i+1]]; j++)
		for (k=0; k<nb[num_lin[i+1]]; k++)
			{
			to_treat[0]=lineage[num_lin[i+1]][j];
			to_treat[1]=lineage[num_lin[i+1]][k];
			if (choices->seqtype) p12=Poisson(to_treat, PROT_POISSON);
			else p12=NC_dist(to_treat, choices->compute_nc);
			if (p12.k<0) pb_k=1;

			for (r=0; r<nb[0]; r++)
				for (rr=0; rr<nb[0]; rr++)
					{
					to_treat[0]=lineage[0][r];
					to_treat[1]=lineage[0][rr];
					if (choices->seqtype) p34=Poisson(to_treat, PROT_POISSON);
					else p34=NC_dist(to_treat, choices->compute_nc);
					if (p34.k<0) pb_k=1;
					weight_pair=weight[0][r]*weight[0][rr]*weight[i+1][j]*weight[i+1][k];

					if (!pb_k)
						{
						if (choices->seqtype) cov=CovPoisson(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34, PROT_POISSON);
						else cov=NC_cov(p12, para[i][r][j], para[i][rr][j], para[i][r][k], para[i][rr][k], p34, choices->compute_nc);
						if (choices->verbose && (j!=k || r!=rr))
							{
							if (cov>=0) fprintf(out_text, "covariance[(%s, %s), (%s, %s)]=%g (sd=%g)\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name, cov, sqrt(cov));
							else fprintf(out_text, "covariance[(%s, %s), (%s, %s)] not computed, saturation\n", lineage[0][r].name, lineage[num_lin[i+1]][j].name, lineage[0][rr].name, lineage[num_lin[i+1]][k].name);
							}
						var_dK+=weight_pair*cov;
						}
					}
			}
	}

Sllcov=0;

for (i=0; i<nb[num_lin[1]]; i++)
	for (j=0; j<nb[num_lin[2]]; j++)
		{
		to_treat[0]=lineage[num_lin[1]][i];
		to_treat[1]=lineage[num_lin[2]][j];
		if (choices->seqtype) p12=Poisson(to_treat, PROT_POISSON);
		else p12=NC_dist(to_treat, choices->compute_nc);
		if (p12.k<0) pb_k=1;

		for (r=0; r<nb[0]; r++)
			for (rr=0; rr<nb[0]; rr++)
				{
				to_treat[0]=lineage[0][r];
				to_treat[1]=lineage[0][rr];
				if (choices->seqtype) p34=Poisson(to_treat, PROT_POISSON);
				else p34=NC_dist(to_treat, choices->compute_nc);
				if (p34.k<0) pb_k=1;

				weight_pair=weight[0][r]*weight[0][rr]*weight[1][i]*weight[2][j];

				if (!pb_k)
					{
					if (choices->seqtype) cov=CovPoisson(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34, PROT_POISSON);
					else cov=NC_cov(p12, para[0][r][i], para[0][rr][i], para[1][r][j], para[1][rr][j], p34, choices->compute_nc);
					Sllcov+=weight_pair*cov;
					}
				}
		}

for (i=0; i<3; i++) weight[i]=De_allocation(weight[i]);
weight=De_allocation(weight);

if (!pb_k)
	{
	if (choices->verbose) fprintf(out_text, "covariance between ingroups=%g (sd=%g)\n", Sllcov, sqrt(Sllcov));
	dK=Kt[0]-Kt[1];
	if (dK==0) pb_k=1;
	}

if (choices->seqtype==NC)
	for (i=0; i<3; i++)
		{
		fprintf(out_text, "mean GC ");
		fprintf(out_text, "%s: ", lineage_name[num_lin[i]]);
		fprintf(out_text, "%.3lf\n", gc_mean[i]);

		if (out_table) fprintf(out_table, "%g\t", gc_mean[i]);
		}
fprintf(out_text, "\n");

if (!pb_k)
	{
	var_dK-=2*Sllcov;
	if (var_dK<=0) pb_k=1;
	}

if (!pb_k)
	{
	fprintf(out_text, "number of compared sites: %.1lf\n", para[0][0][0].l);
	for (i=0; i<2; i++) fprintf(out_text, "mean K %s: %g\n", lineage_name[num_lin[i+1]], Kt[i]);

	fprintf(out_text, "dK: %g\tstandard deviation (sd): %g\tratio dK/sd: %g\n", dK, sqrt(var_dK), Divide(dK, sqrt(var_dK)));
	fprintf(out_text, "exact probability (P): %g\n\n", Exact_probability(Divide(dK, sqrt(var_dK))));

	if (out_table) fprintf(out_table, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t", para[0][0][0].l, Kt[0], Kt[1], dK, sqrt(var_dK), Divide(dK, sqrt(var_dK)), Exact_probability(Divide(dK, sqrt(var_dK))));
	}
else
	{
	if (dK==0)
		{
		fprintf(out_text, "number of compared sites: %.1lf\n", para[0][0][0].l);
		for (i=0; i<2; i++) fprintf(out_text, "mean K %s: %g\n", lineage_name[num_lin[i+1]], Kt[i]);
		fprintf(out_text, "<< no detectable difference, no need to test!\n\n");

		if (out_table) fprintf(out_table, "%g\t%g\t%g\t0\t\t0\t1\t", para[0][0][0].l, Kt[0], Kt[1]);
		}
	else
		{
		fprintf(out_text, "** impossible to compute K. saturation? **\n\n");
		if (out_table) fprintf(out_table, "\t\t\t\t\t\t\t");
		}
	}

for (i=0; i<2; i++)
	{
	for (j=0; j<nb[0]; j++) para[i][j]=De_allocation(para[i][j]);
	para[i]=De_allocation(para[i]);
	}

if (out_table) fprintf(out_table, "\n");
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

struct parameters_lwl Lwl(struct sequences s[2], struct options *choices)
{
static double *tl0[64], *tl1[64], *tl2[64], *tti0[64], *tti1[64], *tti2[64], *ttv0[64], *ttv1[64], *ttv2[64];
int i, j, j3, lgseq, n;
double *rl[21];
char *seq[2];
struct parameters_lwl results;
char *error_message;

if (choices->first)	/*determination of weights of different paths between codons*/
	{
	for (i=0; i<64; i++)
		{
		tl0[i] = Allocation(64, sizeof(double), "tl0[i]");
		tl1[i] = Allocation(64, sizeof(double), "tl1[i]");
		tl2[i] = Allocation(64, sizeof(double), "tl2[i]");
  	 	tti0[i] = Allocation(64, sizeof(double), "tti0[i]");
 		tti1[i] = Allocation(64, sizeof(double), "tti1[i]");
 		tti2[i] = Allocation(64, sizeof(double), "tti2[i]");
 		ttv0[i] = Allocation(64, sizeof(double), "ttv0[i]");
 		ttv1[i] = Allocation(64, sizeof(double), "ttv1[i]");
 		ttv2[i] = Allocation(64, sizeof(double), "ttv2[i]");
		}

	for (i=0; i<21; i++) rl[i] = Allocation(21,sizeof(double), "rl[i]");
	for (i=0; i<21; i++) rl[0][i]=rl[i][0]=0.;

	rl[2][1]=rl[3][1]=rl[3][2]=rl[5][1]=rl[5][3]
          =rl[6][1]=rl[6][3]=rl[6][5]=0.382;
	rl[7][1]=rl[7][3]=rl[7][5]=rl[7][6]=rl[8][1]
          =rl[8][5]=rl[8][6]=rl[8][7]=0.382;
	rl[11][9]=rl[12][9]=rl[13][4]=rl[13][9]=rl[15][4]
           =rl[15][13]=rl[16][14]=rl[16][15]=0.382;
	rl[17][4]=rl[18][4]=rl[18][15]=rl[18][17]=rl[19][4]
           =rl[19][15]=rl[19][16]=rl[20][16]=0.382;
	rl[20][19]=0.382;

	rl[4][1]=rl[4][3]=rl[5][2]=rl[5][4]=rl[6][2]
          =rl[6][4]=rl[7][2]=rl[7][4]=0.343;
	rl[8][2]=rl[8][3]=rl[8][4]=rl[9][4]=rl[9][5]
          =rl[9][6]=rl[9][7]=rl[9][8]=0.343;
	rl[11][4]=rl[11][5]=rl[11][6]=rl[11][7]=rl[11][8]
           =rl[12][4]=rl[12][11]=rl[13][3]=0.343;
	rl[13][5]=rl[13][6]=rl[13][7]=rl[13][8]=rl[13][11]
           =rl[13][12]=rl[14][4]=rl[14][9]=0.343;
	rl[14][11]=rl[14][12]=rl[14][13]=rl[15][3]=rl[15][8]
            =rl[15][9]=rl[15][11]=rl[15][12]=0.343;
	rl[15][14]=rl[16][4]=rl[16][9]=rl[16][12]=rl[16][13]
            =rl[17][3]=rl[17][5]=rl[17][8]=0.343;
	rl[17][13]=rl[17][15]=rl[17][16]=rl[18][1]=rl[18][3]
            =rl[18][5]=rl[18][7]=rl[18][8]=0.343;
	rl[18][13]=rl[18][16]=rl[19][9]=rl[19][12]=rl[19][13]
            =rl[19][14]=rl[19][17]=rl[19][18]=0.343;
	rl[20][4]=rl[20][12]=rl[20][13]=rl[20][14]=rl[20][15]
           =rl[20][18]=0.343;

	rl[4][2]=rl[9][1]=rl[9][2]=rl[9][3]=rl[11][1]
          =rl[11][2]=rl[11][3]=rl[12][5]=0.128;
	rl[12][6]=rl[12][7]=rl[12][8]=rl[13][1]=rl[13][2]
           =rl[13][10]=rl[14][3]=rl[14][5]=0.128;
	rl[14][6]=rl[14][7]=rl[14][8]=rl[14][10]=rl[15][1]
           =rl[15][2]=rl[15][5]=rl[15][6]=0.128;
	rl[15][7]=rl[16][3]=rl[16][5]=rl[16][7]=rl[16][8]
           =rl[16][10]=rl[16][11]=rl[17][1]=0.128;
	rl[17][2]=rl[17][6]=rl[17][7]=rl[17][9]=rl[17][11]
           =rl[17][12]=rl[17][14]=rl[18][2]=0.128;
	rl[18][6]=rl[18][9]=rl[18][11]=rl[18][12]=rl[18][14]
           =rl[19][1]=rl[19][3]=rl[19][5]=0.128;
	rl[19][6]=rl[19][7]=rl[19][8]=rl[19][11]=rl[20][9]
       =rl[20][11]=rl[20][17]=rl[12][3]=0.128;

	rl[10][1]=rl[10][2]=rl[10][3]=rl[10][4]=rl[10][5]
           =rl[10][6]=rl[10][7]=rl[10][8]=0.040;
	rl[10][9]=rl[11][10]=rl[12][1]=rl[12][2]=rl[12][10]
           =rl[14][1]=rl[14][2]=rl[15][10]=0.040;
	rl[16][1]=rl[16][2]=rl[16][6]=rl[17][10]=rl[18][10]
           =rl[19][2]=rl[19][10]=rl[20][1]=0.040;
	rl[20][2]=rl[20][3]=rl[20][5]=rl[20][6]=rl[20][7]
           =rl[20][8]=rl[20][10]=0.040;

	for (i = 1; i <= 20; i++)
		{
  		*(rl[i] + i) = 1.0;
		for (j = i + 1; j <= 20; j++)
		*(rl[i] + j) = *(rl[j] + i);
		}

	PreFastlwl(rl, tl0, tl1, tl2, tti0, tti1, tti2, ttv0, ttv1, ttv2, choices);
	}

for (i=0; i<2; i++) seq[i]=s[i].seq;
lgseq=s[0].lg;

if (Fastlwl(seq, 2, lgseq, tti0, tti1, tti2, ttv0, ttv1, ttv2, tl0, tl1, tl2, &results, choices)<0)
	{
	error_message=Allocation(30+strlen(s[0].name)+strlen(s[1].name), sizeof(char), "error_message");
	sprintf(error_message, "problem with sequences %s and %s\n", s[0].name, s[1].name);
	Leave_program(error_message);
	}

return results;
}

/*****************************************************************************/

/* function of computation of Ks, Ka */ 
/* originally written by Nicolas Galtier */

int Fastlwl(char **seq, int nbseq, int lgseq, double **tti0, double **tti1, double **tti2, double **ttv0, double **ttv1, double **ttv2, double **tl0, double **tl1, double **tl2, struct parameters_lwl *results, struct options *choices)
{

double l[3], a[3], b[3], p[3], q[3],  ti[3], tv[3], cc[3], aaa[3], bb[3], va[3], vb[3],es1,es2,es3,es4,es5,es6,es7,es8;
char ci1, ci2, ci3, cj1, cj2, cj3, cod1[3], cod2[3];
int i, j, ii, jj, nbdiff, cat, pos[3], num1, num2, sat, sat1, sat2;
double *ka[1], *ks[1], *rl[21], *vka[1], *vks[1];

for (i=0; i< nbseq - 1; i++)
	{
	ka[i] = Allocation(2 + 1, sizeof(double), "ka[i]");
	vka[i] = Allocation(2 + 1, sizeof(double), "vka[i]");
	ks[i] = Allocation(2 + 1, sizeof(double), "ks[i]");
	vks[i] = Allocation(2 + 1, sizeof(double), "vks[i]");
	}

sat = sat1 = sat2 = 2;

if ((lgseq%3)!=0)
	{
	fprintf(stderr, "sequence length %d non multiple of 3.\n", lgseq);
	return -1;
	}

for (i = 0; i < nbseq - 1; i++)
	{
	for (j = i + 1; j < nbseq; j++)
		{
		l[0] = l[1] = l[2] = 0;
		ti[0] = ti[1] = ti[2] = tv[0] = tv[1] = tv[2] = 0;

		for (ii = 0; ii < Divide(lgseq, 3); ii++)
			{
			cod1[0] = *(seq[i] + 3 * ii);
			cod1[1] = *(seq[i] + 3 * ii + 1);
			cod1[2] = *(seq[i] + 3 * ii + 2);
			cod2[0] = *(seq[j] + 3 * ii);
			cod2[1] = *(seq[j] + 3 * ii + 1);
			cod2[2] = *(seq[j] + 3 * ii + 2);
			num1 = Num(cod1, choices);
			num2 = Num(cod2, choices);
			l[0] += *(tl0[num1] + num2);
			l[1] += *(tl1[num1] + num2);
			l[2] += *(tl2[num1] + num2);
			ti[0] += *(tti0[num1] + num2);
			ti[1] += *(tti1[num1] + num2);
			ti[2] += *(tti2[num1] + num2);
			tv[0] += *(ttv0[num1] + num2);
			tv[1] += *(ttv1[num1] + num2);
			tv[2] += *(ttv2[num1] + num2);
			}


		for (ii = 0; ii < 3; ii++)
			{
			p[ii] = Divide(ti[ii], l[ii]);
			q[ii] = Divide(tv[ii], l[ii]);
			aaa[ii] = Divide(1., 1 - 2 * p[ii] - q[ii]);
			bb[ii] = Divide(1., 1 - 2 * q[ii]);
			cc[ii] = Divide(aaa[ii] - bb[ii], 2.);	/* Sandrine Hugues did some thorough testing which allowed spotting a bug here: it was written '+' instead of '-'. Changed in version 1.1.4. The final difference is of the order of 1% of the exact probability. */

			if (bb[ii] <= 0)
				{
				b[ii] = 10;
				}
			else
				b[ii] = 0.5 * (double) log(bb[ii]);

			if ((aaa[ii] <= 0) || (bb[ii] <= 0))
				{
				a[ii] = 10;
				}
			else
				a[ii] = 0.5 * (double) log(aaa[ii]) - 0.25 * (double)log(bb[ii]);

			es1=aaa[ii] * aaa[ii] * p[ii] + cc[ii] * cc[ii] * q[ii];
			es2=(aaa[ii] * p[ii] + cc[ii] * q[ii]) * ( aaa[ii] * p[ii] + cc[ii] * q[ii]);

			va[ii] = Divide(aaa[ii]*aaa[ii]*p[ii]+cc[ii]*cc[ii]*q[ii]-(aaa[ii]*p[ii]+cc[ii]*q[ii])*(aaa[ii]*p[ii]+cc[ii]*q[ii]), l[ii]);
			vb[ii] = Divide(bb[ii]*bb[ii]*q[ii]*(1-q[ii]), l[ii]);
			}

		if ((a[1] != 10) && (a[2] != 10) && (b[2] != 10))
			{
			ks[i][j] = Divide(l[1]*a[1]+l[2]*a[2], (l[2]+l[1])) + b[2];
			vks[i][j] = Divide(l[1]*l[1]*va[1]+l[2]*l[2]*va[2], (l[1]+l[2])*(l[1]+l[2])) + vb[2] - Divide(bb[2]*q[2]*(2*aaa[2]*p[2]-cc[2]*(1-q[2])), l[1]+l[2]);
			}
		else
			{
			sat1 = 1;
			vks[i][j]=ks[i][j] = -1;
			}

		if ((a[0] != 10) && (b[0] != 10) && (b[1] != 10))
			{
			ka[i][j] = a[0] + Divide(l[0]*b[0]+l[1]*b[1], l[0]+l[1]);
			vka[i][j] = Divide(l[0]*l[0]*vb[0]+l[1]*l[1]*vb[1],  (l[1]+l[0])*(l[1]+l[0])) + va[0] - Divide(bb[0]*q[0]*(2*aaa[0]*p[0]-cc[0]*(1-q[0])), l[1]+l[0]);
			}
		else
			{
			vka[i][j]=ka[i][j] = -1;
			sat2 = 1;
			}


		results->ks=ks[0][1];
		results->vks=vks[0][1];
		results->ka=ka[0][1];
		results->vka=vka[0][1];
		results->ls=Divide(l[1]*a[1], (a[1]+b[1])) + l[2];
		results->la=Divide(l[1]*b[1], (a[1]+b[1])) + l[0];
		results->a0=a[0];
		results->a2=a[1];
		results->a4=a[2];
		results->b0=b[0];
		results->b2=b[1];
		results->b4=b[2];
		results->vb4=vb[2];
		results->l0=l[0];
		results->l2=l[1];
		results->l4=l[2];

		if (a[2]!=10 && a[1]!=10)
			{
			results->as=Divide(l[2]*b[2]+l[1]*b[1], l[2]+l[1]);
			results->vas=Divide(l[2]*l[2]*vb[2] + l[1]*l[1]*vb[1], (l[1] + l[2])*(l[1]+l[2]));
			}
		else
			{
			results->as=-1;
			results->vas=-1;
			}

		if (b[0]!=10 && b[1]!=10)
			{
			results->ba=Divide(l[0]*b[0]+l[1]*b[1], l[0]+l[1]);
			results->vba=Divide(l[0]*l[0]*vb[0] + l[1]*l[1]*vb[1], (l[1] + l[0])*(l[1]+l[0]));
			}
		else
			{
			results->ba=-1;
			results->vba=-1;
			}
/*
First letter: a=transitions b=transversions k=total
Second letter: s=synonymous, a=non synonymous
- in the function: [0]=non degenerated, [1]=partialy degenerated, [2]=totaly degenerated
- in the structure parameters_lwl: 0=non degenerated, 2=partialy degenerated, 4=totaly degenerated
*/
		}
	}

for (i=0; i< nbseq - 1; i++)
	{
	ka[i] = De_allocation(ka[i]);
	vka[i] = De_allocation(vka[i]);
	ks[i] = De_allocation(ks[i]);
	vks[i] = De_allocation(vks[i]);
	}

if (sat1 == 1)
	sat = 1;
if (sat2 == 1)
	sat = 0;

return sat;
}

/*****************************************************************************/

/* weighting codons */ 
/* originally written by Nicolas Galtier */

int Num(char *cod, struct options *choices)
{
int n1, n2, n3;

n1 = n2 = n3 = 0;

if (cod[0] == 'C')
	n1 = 1;
if (cod[1] == 'C')
	n2 = 1;
if (cod[2] == 'C')
	n3 = 1;
if (cod[0] == 'G')
	n1 = 2;
if (cod[1] == 'G')
	n2 = 2;
if (cod[2] == 'G')
	n3 = 2;
if (cod[0] == 'T')
	n1 = 3;
if (cod[1] == 'T')
	n2 = 3;
if (cod[2] == 'T')
	n3 = 3;

return 16 * n1 + 4 * n2 + n3;
}

/*****************************************************************************/

/* weighting codons */ 
/* originally written by Nicolas Galtier */

int catsite(char c1, char c2, char c3, int i, struct options *choices)
{
/* returns 0 if site i of codon c1c2c3 is non degenerate */
/*         1                              2-fold degenerate */
/*         2                              4-fold degenerate */

if (i == 3)
	{
	if( !choices->code_mt )
		{
		if ( (c1 == 'A') && (c2 == 'T') && (c3 == 'G'))
			return 0;
		if ( (c1 == 'T') && (c2 == 'G') && (c3 == 'A'))
			return 0;
		if ( (c1 == 'T') && (c2 == 'G') && (c3 == 'G'))
			return 0;
		}

	if (c2 == 'C')
		return 2;
	if ((c1 == 'C') && (c2 == 'T'))
		return 2;
	if ((c1 == 'G') && (c2 == 'T'))
		return 2;
	if ((c1 == 'G') && (c2 == 'G'))
		return 2;
	if ((c1 == 'C') && (c2 == 'G'))
		return 2;
	return 1;
	}
else
	if (i == 1)
		{
		if ((c1 == 'C') && (c2 == 'T') && (c3 == 'A'))
			return 1;
		if ((c1 == 'C') && (c2 == 'T') && (c3 == 'G'))
			return 1;
		if ((c1 == 'T') && (c2 == 'T') && (c3 == 'A'))
			return 1;
		if ((c1 == 'T') && (c2 == 'T') && (c3 == 'G'))
			return 1;

		if( !choices->code_mt )
			{
			if ((c1 == 'A') && (c2 == 'G') && (c3 == 'A'))
				return 1;
			if ((c1 == 'A') && (c2 == 'G') && (c3 == 'G'))
				return 1;
			if ((c1 == 'C') && (c2 == 'G') && (c3 == 'A'))
				return 1;
			if ((c1 == 'C') && (c2 == 'G') && (c3 == 'G'))
				return 1;
			}

		return 0;
		}

return 0;
}

/*****************************************************************************/

/* weighting codons */ 
/* originally written by Nicolas Galtier */

char transf(char nt1, char nt2, struct options *choices)
{
if (nt1 == nt2) return 'S';

if ((nt1 == 'A') && (nt2 == 'C'))
	return 'v';
if ((nt1 == 'A') && (nt2 == 'G'))
		return 'i';
if ((nt1 == 'A') && (nt2 == 'T'))
	return 'v';
if ((nt1 == 'G') && (nt2 == 'C'))
	return 'v';
if ((nt1 == 'G') && (nt2 == 'T'))
	return 'v';
if ((nt1 == 'C') && (nt2 == 'T'))
	return 'i';
if ((nt1 == 'C') && (nt2 == 'A'))
	return 'v';
if ((nt1 == 'G') && (nt2 == 'A'))
	return 'i';
if ((nt1 == 'T') && (nt2 == 'A'))
	return 'v';
if ((nt1 == 'C') && (nt2 == 'G'))
	return 'v';
if ((nt1 == 'T') && (nt2 == 'G'))
	return 'v';
if ((nt1 == 'T') && (nt2 == 'C'))
	return 'i';

fprintf(stderr, "Error\n%c, %c\n", nt1, nt2);
return 'E';
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

void titv1(char *cod1, char *cod2, double weight, double *ti, double *tv, double* l, struct options *choices)
{
int i, j, jj;
char a, b, ci1, ci2, ci3, cj1, cj2, cj3;
char transf(char, char, struct options *);
double weight2 = Divide(weight, 2.);

ci1 = cod1[0];
ci2 = cod1[1];
ci3 = cod1[2];
cj1 = cod2[0];
cj2 = cod2[1];
cj3 = cod2[2];

for (i = 0; i <= 2; i++)
	if (cod1[i] != cod2[i])
		{
		l[catsite(ci1, ci2, ci3, i + 1, choices)]+=0.5 * weight;
		l[catsite(cj1, cj2, cj3, i + 1, choices)]+=0.5 * weight;

		a = cod1[i];
		b = cod2[i];
		if (transf(a, b, choices) == 'i')
			{
			ti[catsite(ci1, ci2, ci3, i + 1, choices)] += 0.5 * weight;
			ti[catsite(cj1, cj2, cj3, i + 1, choices)] += 0.5 * weight;
			}
		else
			{
			tv[catsite(ci1, ci2, ci3, i + 1, choices)] += 0.5 * weight;
			tv[catsite(cj1, cj2, cj3, i + 1, choices)] += 0.5 * weight;
			}

		if( choices->code_mt ) continue;  /*there are no problems of TI non-syno and TV syno with code_mt!*/
		
		if (((ci2 == 'T') && (cj2 == 'T')) || ((ci2 == 'G') && (cj2 == 'G')))	/* T or G together in pos 2 of both codons */
			{
/* incrementation removed for readibility */

if (i==0)	/* pos 1 */
	{	
/* all these cases are transitions in a 2-fold non-syno site for the universel code:
they must be removed from the computation of TI 2-fold (ti[1]) and added to computation of TV 2-fold (tv[1])
for code_mt they are non dege sites who have been treated the right way */

	if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'T') && (cj2 == 'G') && (cj3 == 'A'))
		{
		ti[1] -= 0.5 * weight; /* CGA / TGA */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'T') && (cj2 == 'G') && (cj3 == 'G'))
		{
		ti[1] -= 0.5 * weight; /* CGG / TGG */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'G') && (cj2 == 'G') && (cj3 == 'G'))
		{
		ti[1] -= 0.5 * weight; /* AGG / GGG */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'G') && (cj2 == 'G') && (cj3 == 'A'))
		{
		ti[1] -= 0.5 * weight; /* AGA / GGA */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'T') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'A'))
		{
		ti[1] -= 0.5 * weight; /* TGA / CGA */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'T') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'G'))
		{
		ti[1] -= 0.5 * weight; /* TGG / CGG */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'G') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'G'))
		{
		ti[1] -= 0.5 * weight; /* GGG / AGG */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'G') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'A'))
		{
		ti[1] -= 0.5 * weight; /* GGA / AGA */
		tv[1] += 0.5 * weight;
		}


/* all these cases are: 
universel code: TV syno in sites 2-fold they must be removed from computation of TV 2-fold (tv[1]) and added to computation of TI 2-fold (ti[1])
code_mt: TV non syno in site non dege which were correctly counted
*/

	if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'A'))
		{
		tv[1] -= weight; /* CGA / AGA: TV syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'A'))
		{
		tv[1] -= weight; /* AGA / CGA: TV syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'G'))
		{
		tv[1] -= weight; /* CGG / AGG: TV syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'G'))
		{
		tv[1] -= weight; /* AGG / CGG: TV syno code univ, non code mt */
		ti[1] += weight;
		}
	}

if (i==2)	/* pos 3 */
	{

/* all these cases are:
universel code: TV syno in site 2-fold they must be removed from TV 2-fold (iv[1]) and added to TI 2-fold (ti[1])
code_mt: TV non syno in site 2-fold which were normally counted
*/

	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'T'))
		{
		tv[1] -= weight; /* TV ATA / ATT: syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'T') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A'))
		{
		tv[1] -= weight; /* TV ATT / ATA: syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'C'))
		{
		tv[1] -= weight; /* TV ATA / ATC: syno code univ, non code mt */
		ti[1] += weight;
		}
	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'C') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A'))
		{
		tv[1] -= weight; /* TV ATC / ATA: syno code univ, non code mt */
		ti[1] += weight;
		}

/* these 2 cases are:
universel code: TI non syno in site 2-fold they must be removed from TI 2-fold (ti[1]) and added to TV 2-fold (tv[1])
code_mt: TI syno in site 2-fold which were normally counted
*/

	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'G'))
		{
		ti[1] -= 0.5 * weight; /* TI ATA / ATG: non syno code univ, syno code mt */
		tv[1] += 0.5 * weight;
		}
	if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A'))
		{
		ti[1] -= 0.5 * weight; /* TI ATG / ATA: non syno code univ, syno code mt */
		tv[1] += 0.5 * weight;
		}
	}
/* back to normal incrementation */
		}

	}
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

void titv2(char *cod1, char *cod2, double *ti, double *tv, double* l, int *aa, double **rl, int* pos, struct options *choices)
{
	char            codint1[3], codint2[3];
	int             i, j, n, a1, a2, a3, a4, b1, b2, b3,b4;
	double           l1, l2, p1, p2,f1,f2,f3,f4;
	void            titv1(char *, char *, double, double *, double *,double*, struct options *);


	memcpy(codint1, cod1, 3);
	memcpy(codint2, cod1, 3);
	for (i = 0; i < 2; i++)
		{
		if (cod1[i] != cod2[i])
			codint1[i] = cod2[i];
		if (cod1[i] != cod2[i])
			break;
		}
	for (j = i + 1; j <= 2; j++)
		{
		if (cod1[j] != cod2[j])
			codint2[j] = cod2[j];
		if (cod1[j] != cod2[j])
			break;
		}

	
	l1 = *(rl[aa[Num(cod1, choices)]] + aa[Num(codint1, choices)]) * *(rl[aa[Num(codint1, choices)]] + aa[Num(cod2, choices)]);
	l2 = *(rl[aa[Num(cod1, choices)]] + aa[Num(codint2, choices)]) * *(rl[aa[Num(codint2, choices)]] + aa[Num(cod2, choices)]);

	p1 = Divide(l1, l1 + l2);
	p2 = 1 - p1;
	for (i=0;i<3;i++) if (pos[i]==0) n=i+1;
	l[catsite(cod1[0], cod1[1] ,cod1[2], n, choices)]+=0.333333;
	l[catsite(cod2[0], cod2[1] ,cod2[2], n, choices)]+=0.333333;
	l[catsite(codint1[0], codint1[1] ,codint1[2], n, choices)]+=0.333333*p1;
	l[catsite(codint2[0], codint2[1] ,codint2[2], n, choices)]+=0.333333*p2;
	titv1(cod1, codint1, p1, ti, tv,l, choices);
	titv1(cod2, codint1, p1, ti, tv,l, choices);
	titv1(cod1, codint2, p2, ti, tv,l, choices);
	titv1(cod2, codint2, p2, ti, tv,l, choices);


}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

void titv3(char *cod1, char *cod2, double *ti, double *tv, double* l, int *aa, double **rl, struct options *choices)
{

	char           *codint1[6], *codint2[6];
	int             i, j, ii,a,b,c,d,aaa,aab,aac,aad;
	double           like[6], p[6], somli, rlab, rlbc, rlcd;
	void            titv1(char *, char *, double, double *, double *, double*, struct options *);

	for (i = 0; i < 6; i++)
		{
		codint1[i] = Allocation(3, sizeof(char), "codint1[i]");
		codint2[i] = Allocation(3, sizeof(char), "codint2[i]");
		}
	for (i = 0; i < 3; i++)
		{
		for (j = 0; j < 3; j++)
			if (j != i)
				{
				if ((i == 0) || ((i == 1) && (j == 0)))
					{
					ii = 3 * i + j - 1;
					}
				else
					{
					ii = 3 * i + j - 2;
					}
				memcpy(codint1[ii], cod1, 3);
				*(codint1[ii] + i) = cod2[i];
				memcpy(codint2[ii], codint1[ii], 3);
				*(codint2[ii] + j) = cod2[j];
				a=Num(cod1, choices);
				b=Num(codint1[ii], choices);
				c=Num(codint2[ii], choices);
				d=Num(cod2, choices);
				aaa=aa[a];
				aab=aa[b];
				aac=aa[c];
				aad=aa[d];
				rlab=*(rl[aaa]+aab);
				rlbc=*(rl[aab]+aac);
				rlcd=*(rl[aac]+aad);
				like[ii] = rlab*rlbc*rlcd;
				}
		}

	somli = 0;
	for (i = 0; i < 6; i++)
		somli += like[i];
	for (i = 0; i < 6; i++)
		{
		p[i] = Divide(like[i],  somli);
		titv1(cod1, codint1[i], p[i], ti, tv,l, choices);
		titv1(codint1[i], codint2[i], p[i], ti, tv,l, choices);
		titv1(codint2[i], cod2, p[i], ti, tv,l, choices);
		}


}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

void PreFastlwl(double **rl, double **tl0, double **tl1, double **tl2, double **tti0, double **tti1, double **tti2, double **ttv0, double **ttv1, double **ttv2, struct options *choices)
{

	double           l[3], k[3], ti[3], tv[3],cc[3], aaa[3], bb[3], flgseq;
	char             cod1[3], cod2[3];
	int             i, j, ii, jj, nbdiff, cat, pos[3], aa[64], n1, n2, n3;
	void            titv2(char *, char *, double *, double *, double *, int *, double **, int *pos, struct options *choices);
	void            titv3(char *, char *, double *, double *, double *, int *, double **, struct options *choices);
	void            titv1(char *, char *, double, double *, double *, double *, struct options *choices);
	double		 minrl;

/* amino acid code:
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20    0
F W Y H M L I V P  C  A  G  T  S  Q  N  K  R  E  Q  stop
*/

	aa[0] = 17;/* aaa K */
	aa[1] = 16;/* aac N */
	aa[2] = 17;/* aag K */
	aa[3] = 16;/* aat N */
	aa[4] = 13;/* aca T */
	aa[5] = 13;/* acc T */
	aa[6] = 13;/* acg T */
	aa[7] = 13;/* act T */
	if (choices->code_mt)
		aa[8] = 0;/* aga * */
	else
		aa[8] = 18;/* aga R */
	aa[9] = 14;/* agc S */
	if(choices->code_mt)
		aa[10] = 0;/* agg * */
	else
		aa[10] = 18;/* agg R */
	aa[11] = 14;/* agt S */
	if(choices->code_mt)
		aa[12] = 5;/* ata M */
	else
		aa[12] = 7;/* ata I */
	aa[13] = 7;/* atc I */
	aa[14] = 5;/* atg M */
	aa[15] = 7;/* att I */
	aa[16] = 15;
	aa[17] = 4;
	aa[18] = 15;
	aa[19] = 4;
	aa[20] = 9;
	aa[21] = 9;
	aa[22] = 9;
	aa[23] = 9;
	aa[24] = 18;
	aa[25] = 18;
	aa[26] = 18;
	aa[27] = 18;
	aa[28] = 6;
	aa[29] = 6;
	aa[30] = 6;
	aa[31] = 6;
	aa[32] = 19;
	aa[33] = 20;
	aa[34] = 19;
	aa[35] = 20;
	aa[36] = 11;
	aa[37] = 11;
	aa[38] = 11;
	aa[39] = 11;
	aa[40] = 12;
	aa[41] = 12;
	aa[42] = 12;
	aa[43] = 12;
	aa[44] = 8;
	aa[45] = 8;
	aa[46] = 8;
	aa[47] = 8;
	aa[48] = 0;/* taa * */
	aa[49] = 3;/* tac Y */
	aa[50] = 0;/* tag * */
	aa[51] = 3;/* tat Y */
	aa[52] = 14;/* tca S */
	aa[53] = 14;/* tcc S */
	aa[54] = 14;/* tcg S */
	aa[55] = 14;/* tct S */
	if(choices->code_mt)
		aa[56] = 2;/* tga W */
	else
		aa[56] = 0;/* tga * */
	aa[57] = 10;/* tgc */
	aa[58] = 2;/* tgg W */
	aa[59] = 10;/* tgt */
	aa[60] = 6;/* tta */
	aa[61] = 1;/* ttc */
	aa[62] = 6;/* ttg */
	aa[63] = 1;/* ttt */

/* added by Manolo Gouy */
/* compute minrl = minimum value of table rl */
minrl=rl[1][1];
for(i=1; i<=20; i++)
	for(j=i+1; j<=20; j++)
		if(rl[i][j] < minrl ) minrl=rl[i][j];
/* load rl[0][i] and rl[i][0] with minrl for aa = stop */
for(i= 0; i<=20; i++)
	rl[0][i] = rl[i][0] = minrl;


	for (i = 0; i < 64; i++)
		{
		for (j = i; j < 64; j++)
			{
			for(ii=0;ii<3;ii++)
				{
				l[ii]=ti[ii]=tv[ii]=0;
				}

			n1 = Divide(i, 16);
			n2 = Divide(i - 16 * n1, 4);
			n3 = i - 16 * n1 - 4 * n2;
			cod1[0] = 'A';
			if (n1 == 1)
				cod1[0] = 'C';
			if (n1 == 2)
				cod1[0] = 'G';
			if (n1 == 3)
				cod1[0] = 'T';
			cod1[1] = 'A';
			if (n2 == 1)
				cod1[1] = 'C';
			if (n2 == 2)
				cod1[1] = 'G';
			if (n2 == 3)
				cod1[1] = 'T';
			cod1[2] = 'A';
			if (n3 == 1)
				cod1[2] = 'C';
			if (n3 == 2)
				cod1[2] = 'G';
			if (n3 == 3)
				cod1[2] = 'T';

			n1 = Divide(j, 16);
			n2 = Divide(j - 16 * n1, 4);
			n3 = j - 16 * n1 - 4 * n2;
			cod2[0] = 'A';
			if (n1 == 1)
				cod2[0] = 'C';
			if (n1 == 2)
				cod2[0] = 'G';
			if (n1 == 3)
				cod2[0] = 'T';
			cod2[1] = 'A';
			if (n2 == 1)
				cod2[1] = 'C';
			if (n2 == 2)
				cod2[1] = 'G';
			if (n2 == 3)
				cod2[1] = 'T';
			cod2[2] = 'A';
			if (n3 == 1)
				cod2[2] = 'C';
			if (n3 == 2)
				cod2[2] = 'G';
			if (n3 == 3)
				cod2[2] = 'T';

			nbdiff = 0;
			pos[0] = pos[1] = pos[2] = 0;
			if (cod1[0] != cod2[0])
				{
				nbdiff++;
				pos[0] = 1;
				}
			if (cod1[1] != cod2[1])
				{
				nbdiff++;
				pos[1] = 1;
				}
			if (cod1[2] != cod2[2])
				{
				nbdiff++;
				pos[2] = 1;
				}
			if (nbdiff != 2)
				for (jj = 0; jj < 3; jj++)
					if (pos[jj] == 0)
						{
						l[catsite(cod1[0], cod1[1], cod1[2], jj + 1, choices)] += 0.5;
						l[catsite(cod2[0], cod2[1], cod2[2], jj + 1, choices)] += 0.5;
						}
			if (nbdiff == 1)
				titv1(cod1, cod2, 1.0, ti, tv, l, choices);
			if (nbdiff == 2)
				titv2(cod1, cod2, ti, tv, l, aa, rl, pos, choices);
			if (nbdiff == 3)
				titv3(cod1, cod2, ti, tv, l, aa, rl, choices);
			
			*(tl0[i]+j)=*(tl0[j]+i)=l[0];
			*(tl1[i]+j)=*(tl1[j]+i)=l[1];
			*(tl2[i]+j)=*(tl2[j]+i)=l[2];
			*(tti0[i]+j)=*(tti0[j]+i)=ti[0];
			*(tti1[i]+j)=*(tti1[j]+i)=ti[1];
			*(tti2[i]+j)=*(tti2[j]+i)=ti[2];
			*(ttv0[i]+j)=*(ttv0[j]+i)=tv[0];
			*(ttv1[i]+j)=*(ttv1[j]+i)=tv[1];
			*(ttv2[i]+j)=*(ttv2[j]+i)=tv[2];

			}
		}
	
return;	
}

/*****************************************************************************/

/*
                 /\
                /  \
               /    \
             0/      \0'
             /\      /\
            /  \    /  \
           1    2  3    4

cov(K13, K24)=var(K00')
*/

double CovKs(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34)
{
double cov;
double l2, l4, A2, A4, B2, B4, q2, p2, q4, p4, vA2, vA4, vB4;	/*parameters of branch 00' of the tree*/
double a, b, c;	/* computational intermediairies of equations 3 and 4 of Li 1993*/

l2=Divide(p13.l2 + p14.l2 + p23.l2 + p24.l2, 4.);
l4=Divide(p13.l4 + p14.l4 + p23.l4 + p24.l4, 4.);

A2=Divide(p13.a2 + p14.a2 + p23.a2 + p24.a2 - 2*p12.a2 - 2*p34.a2, 4.);
A4=Divide(p13.a4 + p14.a4 + p23.a4 + p24.a4 - 2*p12.a4 - 2*p34.a4, 4.);
B2=Divide(p13.b2 + p14.b2 + p23.b2 + p24.b2 - 2*p12.b2 - 2*p34.b2, 4.);
B4=Divide(p13.b4 + p14.b4 + p23.b4 + p24.b4 - 2*p12.b4 - 2*p34.b4, 4.);

q2=Divide(1. - exp(-2.*B2), 2.);
p2=Divide(1. - exp(-2.*A2 + 0.5*(double)log(1 - 2*q2)) - q2, 2.);

q4=Divide(1. - exp(-2.*B4), 2.);
p4=Divide(1. - exp(-2.*A4 + 0.5*(double)log(1 - 2*q4)) - q4, 2.);

a=Divide(1., 1. - 2.*p2 - q2);
b=Divide(1., 1. - 2.*q2);
c=Divide(a-b, 2.);

vA2=Divide(a*a*p2 + c*c*q2 - (a*p2 +c*q2)*(a*p2 +c*q2), l2);

a=Divide(1., 1. - 2.*p4 - q4);
b=Divide(1., 1. - 2.*q4);
c=Divide(a - b, 2.);

vA4=Divide(a*a*p4 + c*c*q4 - (a*p4 +c*q4)*(a*p4 +c*q4), l4);

vB4=Divide(b*b*q4*(1 - q4), l4);

cov=Divide(l2*l2*vA2 + l4*l4*vA4, (l2 + l4)*(l2 + l4)) + vB4 - Divide(b*q4*(2*a*p4 - c*(1 - q4)), l2 + l4);
/*a, b, c previously computed on positions 4*/

return cov;
}

/*****************************************************************************/

double CovKa(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34)
{
double cov;
double l2, l0, A2, A0, B2, B0, q2, p2, q0, p0, vA0, vB2, vB0;	/*parameters of branch 00' of the tree*/
double a, b, c;	/*computational intermediairies of equations 3 and 4 of Li 1993*/

l2=(p13.l2 + p14.l2 + p23.l2 + p24.l2)/4.;
l0=(p13.l0 + p14.l0 + p23.l0 + p24.l0)/4.;

A2=(p13.a2 + p14.a2 + p23.a2 + p24.a2 - 2*p12.a2 - 2*p34.a2)/4.;
A0=(p13.a0 + p14.a0 + p23.a0 + p24.a0 - 2*p12.a0 - 2*p34.a0)/4.;
B2=(p13.b2 + p14.b2 + p23.b2 + p24.b2 - 2*p12.b2 - 2*p34.b2)/4.;
B0=(p13.b0 + p14.b0 + p23.b0 + p24.b0 - 2*p12.b0 - 2*p34.b0)/4.;

q2=(1. - exp(-2.*B2))/2.;
p2=(1. - exp(-2.*A2 + 0.5*(double)log(1 - 2*q2)) - q2)/2.;

q0=(1. - exp(-2.*B0))/2.;
p0=(1. - exp(-2.*A0 + 0.5*(double)log(1 - 2*q0)) - q0)/2.;

if (((1. - 2.*p2 - q2)==0) || ((1. - 2.*q2)==0)) return -1;

a=Divide(1., 1. - 2.*p2 - q2);
b=Divide(1., 1. - 2.*q2);
c=(a-b)/2.;

vB2=Divide(b*b*q2*(1 - q2), l2);

if (((1. - 2.*p0 - q0)==0) || ((1. - 2.*q0)==0)) return -1;

a=Divide(1., 1. - 2.*p0 - q0);
b=Divide(1., 1. - 2.*q0);
c=(a - b)/2.;

vA0=Divide(a*a*p0 + c*c*q0 - (a*p0 +c*q0)*(a*p0 +c*q0), l0);

vB0=Divide(b*b*q0*(1 - q0), l0);

cov=vA0 + Divide(l0*l0*vB0 + l2*l2*vB2, (l0+l2)*(l0+l2)) - Divide(b*q0*(2*a*p0 - c*(1-q0)), l0+l2);
/*a, b, c previously computed on positions 0*/

return cov;
}

/*****************************************************************************/

double CovAs(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34)
{
double cov;
double l2, l4, A2, A4, B2, B4, q2, p2, q4, p4, vA2, vA4;	/*parameters of branch 00' of the tree*/
double a, b, c;	/*computational intermediairies of equations 3 and 4 of Li 1993*/

l2=(p13.l2 + p14.l2 + p23.l2 + p24.l2)/4.;
l4=(p13.l4 + p14.l4 + p23.l4 + p24.l4)/4.;

A2=(p13.a2 + p14.a2 + p23.a2 + p24.a2 - 2*p12.a2 - 2*p34.a2)/4.;
A4=(p13.a4 + p14.a4 + p23.a4 + p24.a4 - 2*p12.a4 - 2*p34.a4)/4.;
B2=(p13.b2 + p14.b2 + p23.b2 + p24.b2 - 2*p12.b2 - 2*p34.b2)/4.;
B4=(p13.b4 + p14.b4 + p23.b4 + p24.b4 - 2*p12.b4 - 2*p34.b4)/4.;

q2=(1. - exp(-2.*B2))/2.;
p2=(1. - exp(-2.*A2 + 0.5*(double)log(1 - 2*q2)) - q2)/2.;

q4=(1. - exp(-2.*B4))/2.;
p4=(1. - exp(-2.*A4 + 0.5*(double)log(1 - 2*q4)) - q4)/2.;

a=Divide(1., 1. - 2.*p2 - q2);
b=Divide(1., 1. - 2.*q2);
c=(a-b)/2.;

vA2=Divide(a*a*p2 + c*c*q2 - (a*p2 +c*q2)*(a*p2 +c*q2), l2);

a=Divide(1., 1. - 2.*p4 - q4);
b=Divide(1., 1. - 2.*q4);
c=(a - b)/2.;

vA4=Divide(a*a*p4 + c*c*q4 - (a*p4 +c*q4)*(a*p4 +c*q4), l4);

cov=Divide(l2*l2*vA2 + l4*l4*vA4, (l2 + l4)*(l2 + l4)) - Divide(b*q4*(2*a*p4 - c*(1 - q4)), l2 + l4);
/*a, b, c previously computed on positions 4*/

return cov;
}

/*****************************************************************************/

double CovBa(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34)
{
double cov;
double l2, l0, A2, A0, B2, B0, q2, p2, q0, p0, vB2, vB0;	/*parameters of branch 00' of the tree*/
double a, b, c;	/*computational intermediairies of equations 3 and 4 of Li 1993*/

l2=(p13.l2 + p14.l2 + p23.l2 + p24.l2)/4.;
l0=(p13.l0 + p14.l0 + p23.l0 + p24.l0)/4.;

A2=(p13.a2 + p14.a2 + p23.a2 + p24.a2 - 2*p12.a2 - 2*p34.a2)/4.;
A0=(p13.a0 + p14.a0 + p23.a0 + p24.a0 - 2*p12.a0 - 2*p34.a0)/4.;
B2=(p13.b2 + p14.b2 + p23.b2 + p24.b2 - 2*p12.b2 - 2*p34.b2)/4.;
B0=(p13.b0 + p14.b0 + p23.b0 + p24.b0 - 2*p12.b0 - 2*p34.b0)/4.;

q2=(1. - exp(-2.*B2))/2.;
p2=(1. - exp(-2.*A2 + 0.5*(double)log(1 - 2*q2)) - q2)/2.;

q0=(1. - exp(-2.*B0))/2.;
p0=(1. - exp(-2.*A0 + 0.5*(double)log(1 - 2*q0)) - q0)/2.;

a=Divide(1., 1. - 2.*p2 - q2);
b=Divide(1., 1. - 2.*q2);
c=(a-b)/2.;

vB2=Divide(b*b*q2*(1 - q2), l2);

a=Divide(1., 1. - 2.*p0 - q0);
b=Divide(1., 1. - 2.*q0);
c=(a - b)/2.;

vB0=Divide(b*b*q0*(1 - q0), l0);

cov=Divide(l0*l0*vB0 + l2*l2*vB2, (l0+l2)*(l0+l2)) - Divide(b*q0*(2*a*p0 - c*(1-q0)), l0+l2);
/*a, b, c previously computed on positions 0*/

return cov;
}

/*****************************************************************************/

double CovB4(struct parameters_lwl p12, struct parameters_lwl p13, struct parameters_lwl p14, struct parameters_lwl p23, struct parameters_lwl p24, struct parameters_lwl p34)
{
double cov;
double l4, B4, q4;	/*parameters of branch 00' of the tree*/
double b;	/*computational intermediairies of equations 3 and 4 of Li 1993*/

l4=(p13.l4 + p14.l4 + p23.l4 + p24.l4)/4.;

B4=(p13.b4 + p14.b4 + p23.b4 + p24.b4 - 2*p12.b4 - 2*p34.b4)/4.;

q4=(1. - exp(-2.*B4))/2.;

b=Divide(1., 1. - 2.*q4);

cov=Divide(b*b*q4*(1 - q4), l4);

return cov;
}

/*****************************************************************************/

double CovK2(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34)
{
double cov;
double l, A, B, q, p;	/*parameters of branch 00' of the tree*/
double a, b;	/*computational intermediairies*/

l=(p13.l + p14.l + p23.l + p24.l)/4.;

A=(p13.a + p14.a + p23.a + p24.a - 2*p12.a - 2*p34.a)/4.;
B=(p13.b + p14.b + p23.b + p24.b - 2*p12.b - 2*p34.b)/4.;

q=(1. - exp(-2.*B))/2.;
p=(1. - exp(-2.*A + 0.5*(double)log(1 - 2*q)) - q)/2.;

if ((1.-2.*q)<=0 || (1.-2.*p-q)<=0) return -1;

a=Divide(1., 1.-2.*p-q);
b=( a + 1./(1.-2.*q) )/2.;

cov=Divide(a*a*p + b*b*q - (a*p+b*q)*(a*p+b*q), l);

return cov;
}

/*****************************************************************************/

int Count_otu(char *tree)
{
int nb=1;
int i=0;

while (tree[i]!=0)
	{
	if (tree[i]=='(') nb++;
	i++;
	}

return nb;
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */
/* ChaineTOTable. Read c_tree input (string) and write t_tree tree (int**), */
/* branch length lgbi (internal) and lgbp (terminal), bootstrap values, */
/* species names (name) and rooted/unrooted (root-> r (rooted) or n (not)). */

int retder(int *list)
{
int i=0, j;

while (list[i] != 0) i++;

j = *(list + i - 1);
*(list + i - 1) = 0;

return j;
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

void aj(int *list, int nb)
{
int  i=0;

while (list[i] != 0) i++;
*(list + i) = nb;

return;
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */

int Tree_ctot(char *input, int **tree, double *lgbi, double *lgbp, double* bootstrap, char **name, int *root)
{
int   i=0, j, k, fin=0, nbpo=0, nbpf = 0, cptv1 = 0, nbotu, br_ouverte, *listcour, otu = -1, pomoinspf, cpttree=0, t=1;
char  c, cc, cas, dejalu = '\0';
double  f;

sscanf(input+(cpttree++), "%c", &c);

if (c == '[')
	{
	while ((c != ']') && c) sscanf(input+(cpttree++), "%c", &c);
	if (c != ']')
		{
		fprintf(stderr, "Unmatched '[' ']'\n");
		return -1;
		}
  	} 
else
	if (c == '(') cpttree=0;
	else
		{
		fprintf(stderr, "Tree 1st character must be '(' or '['\n and is %c\n", c);
		return -1;
		}

while ((c != ';') && c)
	{
	sscanf(input+(cpttree++), "%c", &c);

	if (c == '(') nbpo++;
	if (c == ')') nbpf++;
	if ((nbpo == nbpf + 1) && (c == ',')) cptv1++;
	}

if (c != ';')
	{
	fprintf(stderr, "';' missing at end of tree\n");
	return -1;
	}

if (nbpo != nbpf)
	{
	fprintf(stderr, "Unmatched parenthesis\n");
	return -1;
	}

if (cptv1 == 1) cas = 'c';
if (cptv1 == 2) cas = 'a';

if ((cptv1!=1) && (cptv1!=2))
	{
	fprintf(stderr, "Bad number of ',' in tree\n");
	return -1;
	}

nbotu = nbpo + 2;

if (root)
	if (cas == 'a') *root=0;
	else *root=1;


if(tree==NULL && lgbi==NULL && lgbp==NULL && bootstrap==NULL && name==NULL) goto end;

if (cas=='c')
	{
	if (lgbp) lgbp[0] = 0.;
	if (name) sprintf(name[0],"ROOT");
	if (tree)
		for(i=0;i<nbotu-3;i++) tree[0][i]=0;

	otu++;
	}

listcour=Allocation(nbotu, sizeof(int), "listcour");

cpttree=0;
sscanf(input+(cpttree++), "%c", &c);

if (c == '[')
	{
	while (c != ']') sscanf(input+(cpttree++), "%c", &c);
	while((c==']') || (c==' ') || (c=='\n') || (c=='\t')) sscanf(input+(cpttree++), "%c", &c);

	if (c!='(') return -1;
	}
else
	while(c!='(') sscanf(input+(cpttree++), "%c", &c);

pomoinspf=1;
for (i = 0; i < nbotu; i++) listcour[i] = 0;
br_ouverte = 0;

for (k = 0; t==1; k++)
	{
	if (dejalu == '\0') sscanf(input+(cpttree++), "%c", &c);
	else
		{
		c = dejalu;
		dejalu = '\0';
		}

	switch (c)
		{
		case ';':
			fin = 1;
			break;

		case ',':
		case '\n':
		case '\t':
		case '\'':
		case ' ':
			break;

		case '(':
			pomoinspf ++;
			br_ouverte++;
			aj(listcour, br_ouverte);
			break;

		case ')':
			pomoinspf--;
			sscanf(input+(cpttree++), "%c", &cc);

			if (cc == ';' || pomoinspf==0)
				{
				fin = 1;
				break;
				}

			j = retder(listcour);

			while (cc=='\n' || cc==' ' || cc=='\t')
				sscanf(input+(cpttree++),"%c",&cc);

			if (strpbrk(input+cpttree-1, "0123456789.")==input+cpttree-1)
				{
				if (bootstrap)
					sscanf(input+cpttree-1, "%le", bootstrap+j-1);
				cpttree+=strspn(input+cpttree-1,".0123456789");
				cc=*(input+cpttree-1);
				while (cc=='\n' || cc==' ' || cc=='\t')
					sscanf(input+(cpttree++),"%c",&cc);
				}

			if (cc == ':')
				{
				while(input[cpttree]==' ') cpttree++;
				sscanf(input+cpttree, "%le", &f);
				cpttree+=strspn(input+cpttree,"-0123456789.");
				if (lgbi) lgbi[j - 1] = f;
				}
			else
				dejalu=cc;

			break;

		default:
			otu++;
			cc = c;
			i = 0;

			while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n') && (cc != ' '))
				{
				if (name && cc != '\'')
					{
					name[otu][i] = cc;
					i++;
					name[otu][i]='\0';
					}
				sscanf(input+(cpttree++), "%c", &cc);
				}

			while(input[cpttree-1]==' ') cpttree++;

			cc=input[cpttree-1];

			if (cc == ':')
				{
				while(input[cpttree]==' ') cpttree++; 
				sscanf(input+(cpttree), "%le", &f);
				cpttree+=strspn(input+cpttree,"-0123456789.e");
				if (lgbp) lgbp[otu] = f;
				} 
			else
				dejalu = cc;

			for (i = 0; i < nbotu - 3; i++)
				if (tree) tree[otu][i] = 0;

			for (i = 0; i < nbotu; i++)
				if (tree && (listcour[i]!=0))
					tree[otu][listcour[i] - 1] = 1;
			}

		if (fin == 1)
			break;
		}

listcour=De_allocation(listcour);

end:
	if (cas=='a') return nbotu;
	else return (nbotu-1);
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */
/* Arbitrarily roots tree. */

void Ctree_root(char* ctree)
{
int i=0, place, diffp=0, flag=0;

while(ctree[i])
	{
	if (ctree[i]=='(') diffp++;
	if (ctree[i]==')') diffp--;
	if (diffp==1 && ctree[i]==',') flag++;
	if (flag==2)
		{
		place=i;
		break;
		}
	i++;
	}

for (i=(int)strlen(ctree); i>=place; i--) ctree[i+2]=ctree[i];

ctree[place+1]=')';

for (i=place-1; i>=0; i--) ctree[i+1]=ctree[i];
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */
/* Return 1 if c_tree carbre is rooted, 0 if unrooted, -1 if problem */

int Rooted(char* carbre)
{
int i=0, cpt=0, cptv=0;

while((carbre[i]==' ' || carbre[i]=='\n' || carbre[i]=='\t') && carbre[i]) i++;

if(carbre[i]!='[' && carbre[i]!='(')
	{
	printf("Tree first char must be ( or [\n"); 
	return -1; 
	}

if(carbre[i]=='[')
	while(carbre[i]!=']' && carbre[i]) i++;

if(!carbre[i])
	{
	printf("Unmatched '[' ']'\n");
	return -1;
	}

while(carbre[i]!='(' && carbre[i]) i++;

if(!carbre[i])
	{
	fprintf(stderr, "No initial parenthesis\n");
	return -1;
	}

while(carbre[i]!=';' && carbre[i])
	{
	if (carbre[i]=='(') cpt++;
	if (carbre[i]==')') cpt--;
	if (carbre[i]==',' && cpt==1) cptv++;
	i++;
	}
 
if (cptv==2) return 0;
if (cptv==1) return 1;
return -1;
}

/*****************************************************************************/

double **Weighting(char *tree_parenthesis, struct sequences **lineage, int *nb, int *num_lin)
{
char **names_tree, **names_partial;
int **tab_tree, nb_otu, nb_lineages, **no_root, **tab_partial, nb_partial;
double **weights;
int i, j, k, l, m, n;
int *ok_tree, **ok_seq;
int num_outg=-1, root;
int otu_removed=0, n0, n1, column_to_remove;
double sum;
int lig;
int *outgroup;
int *weight_node, nb_node;	/*compute number of branches per node*/
char *error_message;

root=Rooted(tree_parenthesis);
if (!root) Ctree_root(tree_parenthesis);

nb_otu=Count_otu(tree_parenthesis);
nb_lineages=nb[0]+nb[num_lin[1]]+nb[num_lin[2]];

names_tree=Allocation(nb_otu+1, sizeof(char*), "names_tree");
tab_tree=Allocation(nb_otu+1, sizeof(int*), "tab_tree");
no_root=Allocation(nb_otu, sizeof(int*), "no_root");

for (i=0; i<nb_otu+1; i++)
	{
	names_tree[i]=Allocation(50, sizeof(char), "names_tree[i]");
	tab_tree[i]=Allocation(nb_otu-2, sizeof(int), "tab_tree[i]");
	}

for (i=0; i<nb_otu; i++)
	no_root[i]=Allocation(nb_otu-3, sizeof(int), "no_root[i]");

outgroup=Allocation(nb_otu, sizeof(int), "outgroup");

Tree_ctot(tree_parenthesis, tab_tree, NULL, NULL, NULL, names_tree, &root);

Polytomies(tree_parenthesis, tab_tree, names_tree, nb_otu);

ok_tree=Allocation(nb_otu+1, sizeof(int), "ok_tree");
ok_tree[0]=1;
ok_seq=Allocation(3, sizeof(int*), "ok_seq");

for (i=0; i<3; i++)	/*3 lineages*/
	{
	ok_seq[i]=Allocation(nb[num_lin[i]], sizeof(int), "ok_seq[i]");
	for (j=0; j<nb[num_lin[i]]; j++)	/*number of species per lineage*/
		for (k=1; k<nb_otu+1; k++)	/*species of the tree*/
			if (strcmp(names_tree[k], lineage[num_lin[i]][j].name)==0)
				{
				ok_seq[i][j]=1;
				ok_tree[k]=i+1;
				break;
				}
	}

for (i=0; i<3; i++)	/*3 lineages*/
	for (j=0; j<nb[num_lin[i]]; j++)	/*number of species per lineage*/
		if (!ok_seq[i][j])
			{
			error_message=Allocation(40+strlen(lineage[num_lin[i]][j].name), sizeof(char), "error_message");
			sprintf(error_message, "Sequence %s is absent from the tree", lineage[num_lin[i]][j].name);
			Leave_program(error_message);
			}

for (i=1; i<nb_otu+1; i++)
	{
	for (j=0; j<nb[num_lin[0]]; j++)
		if (strcmp(names_tree[i], lineage[num_lin[0]][j].name)==0)
			{
			outgroup[i-1]=1;
			if (num_outg==-1) num_outg=i;
			else num_outg=-2;		/*num_outg is -2 if there are several sequences in the outgroup*/
			break;
			}
	}

/*change in root position, between outgroup and ingroups*/
if (num_outg>=0)
	{
	if (root)
		{
		Unroot(tab_tree, no_root, nb_otu);
		num_outg--;
		}
	Root(no_root, tab_tree, num_outg, nb_otu);
	}
else
	{
	if (root) Unroot(tab_tree, no_root, nb_otu);

	/*find the internal branch separating the outgroups from the rest*/
	for (i=0; i<nb_otu-3; i++)	/*go through all the internal branches*/
		{
		num_outg=i;
		k=-1;
		for (j=0; j<nb_otu; j++)	/*all taxa*/
			if (ok_tree[j+1])	/*only use taxa which will remain in the analysis, j+1 because the root has been removed*/
				{
				if (outgroup[j])
					{
					if (k<0)
						k=no_root[j][i];
					else
						if (no_root[j][i]!=k)
							{
							num_outg=-1;
							break;
							}
					}
				else
					{
					if (k<0)
						k=1-no_root[j][i];
					else
						if (no_root[j][i]==k)
							{
							num_outg=-1;
							break;
							}
					}
				}
		if (num_outg>=0) break;
		}

	if (num_outg<0) Leave_program("No internal branch separating the outgroups from other sequences found.");

	Root_bi(no_root, tab_tree, num_outg, nb_otu);
	}

/*removing excluded taxa*/

for (i=0; i<nb_otu+1; i++)
	if (!ok_tree[i])
		{
		Remove_line(i, nb_otu, tab_tree, names_tree);
		Remove_column(nb_otu, tab_tree);

		for (j=i; j<nb_otu; j++) ok_tree[j]=ok_tree[j+1];
		ok_tree[nb_otu]=0;

		nb_otu--;
		i--;
		}

weights=Allocation(3, sizeof(double*), "weights");

for (i=0; i<3; i++)
	{
	tab_partial=Allocation(nb_otu+1, sizeof(int*), "tab_partial");
	for (j=0; j<nb_otu+1; j++) tab_partial[j]=Allocation(nb_otu-2, sizeof(int), "tab_partial[i]");
	names_partial=Allocation(nb_otu+1, sizeof(char*), "names_partial");
	for (j=0; j<nb_otu+1; j++) names_partial[j]=Allocation(50, sizeof(char), "names_partial[i]");

	for (j=0; j<nb_otu+1; j++)
		{
		for (k=0; k<nb_otu-2; k++) tab_partial[j][k]=tab_tree[j][k];
		for (k=0; k<(int)strlen(names_tree[j])+1; k++) names_partial[j][k]=names_tree[j][k];
		}
	nb_partial=nb_otu;

	if (i)
		{
		otu_removed=0;

		for (j=1; j<nb_partial+1; j++)
			if ((ok_tree[j+otu_removed]-1)==(3-i))	/*3-i: 1->2 et 2->1*/
				{
				Remove_line(j, nb_partial, tab_partial, names_partial);
				Remove_column(nb_partial, tab_partial);

				nb_partial--;
				otu_removed++;
				j--;
				}
		}

	if (nb_partial-2) weight_node=Allocation(nb_partial-2, sizeof(int), "weight_node");
	else weight_node=NULL;

	for (j=0; j<nb_partial-2; j++)
		for (k=1; k<nb_partial+1; k++)
			weight_node[j]+=tab_partial[k][j];

/*sort nodes by number of branches*/
	for (j=0; j<nb_partial-3; j++)
		for (k=j+1; k<nb_partial-2; k++)
			if (weight_node[k]>weight_node[j])
				{
				m=weight_node[k];
				weight_node[k]=weight_node[j];
				weight_node[j]=m;

				for (l=1; l<nb_partial+1; l++)
					{
					m=tab_partial[l][k];
					tab_partial[l][k]=tab_partial[l][j];
					tab_partial[l][j]=m;
					}
				}

	for (j=0; j<nb_partial-2; j++)	/*all internal branches*/
		if (weight_node[j])
			{
			for (l=j+1; l<nb_partial-2; l++)	/*all branches which may be included in the j-th*/
				{
				nb_node=1;

				for (k=1; k<nb_partial+1; k++)	/*all taxa*/
					if (tab_partial[k][j]==1 && tab_partial[k][l]==1)
						{
						for (n=0; n<nb_partial-2; n++) tab_partial[k][n]*=2;

						for (m=k+1; m<nb_partial+1; m++)
							if (tab_partial[m][l]==1)
								{
								nb_node++;
								for (n=0; n<nb_partial-2; n++) tab_partial[m][n]*=2;	/*so that subdivisions of this lineage should not be used*/
								}
						break;
						}
				weight_node[j]-=(nb_node-1);
				}

			for (l=j+1; l<nb_partial-2; l++)
				for (k=1; k<nb_partial+1; k++)
					if (tab_partial[k][l]>1) tab_partial[k][l]=1;
			}

	weights[i]=Allocation(nb[num_lin[i]], sizeof(double), "weights[i]");
	sum=0;

	for (j=0; j<nb[num_lin[i]]; j++)
		{
		weights[i][j]=1.;

		for (k=1; k<nb_partial+1; k++)
			if (strcmp(lineage[num_lin[i]][j].name, names_partial[k])==0)
				{
				for (l=0; l<nb_partial-2; l++)
					if (tab_partial[k][l])
						weights[i][j]*=Divide(1., weight_node[l]);
				break;
				}
		sum+=weights[i][j];
		}

	for (j=0; j<nb[num_lin[i]]; j++) weights[i][j]*=Divide(1., sum);

	for (j=0; j<nb_otu+1; j++) tab_partial[j]=De_allocation(tab_partial[j]);
	tab_partial=De_allocation(tab_partial);
	for (j=0; j<nb_otu+1; j++) names_partial[j]=De_allocation(names_partial[j]);
	names_partial=De_allocation(names_partial);
	weight_node=De_allocation(weight_node);
	}

for (i=0; i<nb_otu+1; i++)
	{
	names_tree[i]=De_allocation(names_tree[i]);
	tab_tree[i]=De_allocation(tab_tree[i]);
	}
for (i=0; i<nb_otu; i++) no_root[i]=De_allocation(no_root[i]);
names_tree=De_allocation(names_tree);
tab_tree=De_allocation(tab_tree);
no_root=De_allocation(no_root);
outgroup=De_allocation(outgroup);
ok_tree=De_allocation(ok_tree);
for (i=0; i<3; i++) ok_seq[i]=De_allocation(ok_seq[i]);
ok_seq=De_allocation(ok_seq);

return weights;
}

/*****************************************************************************/

void Remove_line(int remove, int nb_otu, int **tab_tree, char **names_tree)
{
int i, j;

for (i=remove; i<nb_otu; i++)
	{
	for (j=0; j<nb_otu-2; j++) tab_tree[i][j]=tab_tree[i+1][j];
	strcpy(names_tree[i], names_tree[i+1]);
	}

for (j=0; j<nb_otu-2; j++) tab_tree[nb_otu][j]=-1;

names_tree[nb_otu]=NULL;
}

/*****************************************************************************/

void Remove_column(int nb_otu, int **tab_tree)
{ 
int i, j, k;
int column_to_remove=-1, *n0, *n1;

n0=Allocation(nb_otu-2, sizeof(int), "n0");
n1=Allocation(nb_otu-2, sizeof(int), "n1");

/*search a non informative column (only one 1 or only one 0)*/
for (i=0; i<nb_otu-2; i++)	/*columns*/
	{
	for (j=0; j<nb_otu; j++)	/*lines*/
		{
		if (tab_tree[j][i]) n1[i]++;
		else n0[i]++;

		if ((n0[i]>1) && (n1[i]>1)) break;
		}

	if ((n0[i]==1) || (n1[i]==1))
		{
		column_to_remove=i;
		break;
		}
	}

if (column_to_remove<0)
/*searche for a column present twice*/
	for (i=0; i<nb_otu-2; i++)	/*columns*/
		if (n1[i])	/*except 0 only columns (due to polytomies)*/
			{
			for (j=i+1; j<nb_otu-2; j++)	/*the other columns, to compare*/
				{
				column_to_remove=i;
				for (k=0; k<nb_otu; k++)	/*lines*/
					if (tab_tree[k][i]!=tab_tree[k][j])
						{
						column_to_remove=-1;
						break;
						}
				if (column_to_remove==i) break;
				}
			if (column_to_remove==i) break;
			}

if (column_to_remove<0)
/*searche for a column with only 0*/
	for (i=0; i<nb_otu-2; i++)	/*columns*/
		if (!n1[i])
			{
			column_to_remove=i;
			break;
			}

if (column_to_remove<0) Leave_program("problem removing columns.");

for (i=0; i<nb_otu; i++)	/*lines*/
	{
	for (j=column_to_remove; j<nb_otu-2-1; j++)	/*columns*/
		tab_tree[i][j]=tab_tree[i][j+1];
	tab_tree[i][nb_otu-3]=-1;
	}

n0=De_allocation(n0);
n1=De_allocation(n1);
}

/*****************************************************************************/

void Polytomies(char *text, int **tab, char **names, int nb_otu)
{
int r1=1, r2, i, j, k;
int nb_parenthesis, nb_name, *tab_polyt, first_polyt, first_non_polyt;
char **name_polyt;

tab_polyt=Allocation(nb_otu+1, sizeof(int), "tab_polyt");
name_polyt=Allocation(nb_otu+1, sizeof(char*), "name_polyt");
for (i=0; i<nb_otu+1; i++) name_polyt[i]=Allocation(50, sizeof(char), "name_polyt[i]");

while (text[r1+5] && text[r1+5]!=';')
	if (text[r1-1]==':' && Is_text_zero(text, r1))
		{
		r2=r1-1;
		while (r2>=0 && text[r2]!=')') r2--;

		nb_parenthesis=1;
		while (r2>=0 && nb_parenthesis)
			{
			r2--;
			if (text[r2]==')') nb_parenthesis++;
			if (text[r2]=='(') nb_parenthesis--;
			}

		/*r2 and r1 coordinates between which the polytomy is*/
		/*find names in this interval*/

		nb_name=0;

		for (i=r2; i<r1; i++)
			{
			j=0;
			if ((text[i]=='(' || text[i]==',') && text[i+1]!='(') j=i+1;
			while (text[j]==' ') j++;
			if (text[j]=='(') j=0;

			if (j)
				{
				k=0;
				while (text[j]!=':')
					{
					name_polyt[nb_name][k]=text[j];
					k++;
					j++;
					}
				name_polyt[nb_name][k]=0;
				nb_name++;
				}
			}

		/*find and put to 0 the column which distinguishes these names and them alone from the rest*/

		first_polyt=-1;
		first_non_polyt=-1;
		for (i=1; i<nb_otu+1; i++)
			{
			tab_polyt[i]=0;
			for (j=0; j<nb_name; j++)
				if (strcmp(name_polyt[j], names[i])==0)
					{
					tab_polyt[i]=1;
					if (first_polyt<0) first_polyt=i;
					break;
					}
			if (first_non_polyt<0 && tab_polyt[i]==0) first_non_polyt=i;
			}

		for (i=0; i<nb_otu-2; i++)
			{
			for (j=1; j<nb_otu+1; j++)
				if ((tab_polyt[j] && tab[j][i]!=tab[first_polyt][i]) || (!tab_polyt[j] && tab[j][i]!=tab[first_non_polyt][i]))
					break;

			if (j==(nb_otu+1)) for (k=0; k<nb_otu+1; k++) tab[k][i]=0;
			}

		r1++;
		}
	else r1++;

tab_polyt=De_allocation(tab_polyt);
for (i=0; i<nb_otu+1; i++) name_polyt[i]=De_allocation(name_polyt[i]);
name_polyt=De_allocation(name_polyt);
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */
/*
The function "root" puts in the "r-tree" table the topology corresponding to that of
the tree once rooted on the peripherical branch of taxon "otu".
"notu" must contain the size of the tree = the number of lines = the number of columns+3.
The allocation of "r-tree" must be done at input.
The order of taxa is not modified.
BEWARE: the list of names does not contain the now necessary "ROOT".
*/

void Root(int** tree, int** r_tree, int otu, int notu)
{
int i,j;

for(i=0;i<notu-2;i++) r_tree[0][i]=0;
	
for(i=0;i<notu-3;i++)
	{
	if (tree[otu][i]==1)
		{
		for(j=1;j<notu+1;j++)
			{
			if(tree[j-1][i]==0) r_tree[j][i]=1;
			else r_tree[j][i]=0;
			}
		}
	else
		{
		for(j=1;j<notu+1;j++) r_tree[j][i]=tree[j-1][i];
		}
	}
	
for(i=0;i<notu+1;i++) r_tree[i][notu-3]=1;

r_tree[0][notu-3]=r_tree[otu+1][notu-3]=0;

return;
}

/*****************************************************************************/

/* originally written by Nicolas Galtier */
/*
The function "root_bi" puts in the table "r_tree" the topology of the tree once rooted
on internal branch number "bi".
"notu" must contain the size of the tree = the number of lines = the number of columns+3.
The allocation of "r_tree" must be done at input.
*/

void Root_bi(int** tree, int** r_tree, int bi, int notu)
{
int i, j, flag10, flag01, sumi, sumbi=0;

for(i=0;i<notu;i++) sumbi+=tree[i][bi];

for(i=0;i<notu-2;i++) r_tree[0][i]=0;

for(i=0;i<notu-3;i++)
	{
	if (i!=bi)
		{
		sumi=0;

		for(j=0;j<notu;j++) sumi+=tree[j][i];

		flag10=flag01=0;

		for(j=0;j<notu;j++)
			{
			if(tree[j][i]==0 && tree[j][bi]==1)
				{ 
				flag01++;
				break;
				}
			}

		for(j=0;j<notu;j++)
			{
			if(tree[j][i]==1 && tree[j][bi]==0)
				{ 
				flag10++;
				break;
				}
			}

		if ((flag01==1 && flag10==1 && notu-sumi>sumbi) || (!(flag01==1 && flag10==1) && sumi<sumbi) )
			{
			for(j=0;j<notu;j++)
				{
				r_tree[j+1][i]=tree[j][i];
				}
			}
		else
			{
			for(j=0;j<notu;j++)
				{
				if(tree[j][i]==0) r_tree[j+1][i]=1;
				if(tree[j][i]==1) r_tree[j+1][i]=0;
				}
			}

		}
	else
		{
		for(j=0;j<notu;j++)
			r_tree[j+1][i]=tree[j][i];
		}
	}

for(j=1;j<notu+1;j++)
	{
	if (tree[j-1][bi]==1) r_tree[j][notu-3]=0;
	else r_tree[j][notu-3]=1;
	}
}
	
/*****************************************************************************/

/* originally written by Nicolas Galtier */
/*
The function "unroot" unroots the rooted tree "treer" and writes the resulting
unrooted tree in "treenr".
"notu" is the number of taxa = the number of lines of "treenr" = the number of
lines of "treer"-1.
The allocation of "treenr" must be done at input.
The order of taxa is not changed.
BEWARE: the list of names still contains the now unwanted "ROOT".
*/

int Unroot(int** treer, int** treenr, int notu)
{
int i, j, j1, j2, br_to_remove=-1, sum1, sum2, flag;

for(j=0;j<notu-2;j++)
	if (treer[0][j]!=0) return 1;
	
/* determination of the internal branch which must disapear */

for(j1=0;j1<notu-3;j1++)
	{
	sum1=0;
	for(i=0;i<notu+1;i++) sum1 += treer[i][j1];

	if (sum1==notu-1)
		{
		br_to_remove=j1;
		break;
		}

	for(j2=j1+1;j2<notu-2;j2++)
		{
		sum2=0;
		for(i=0;i<notu+1;i++) sum2+=treer[i][j2];

		if (sum1==notu-sum2)
			{
			flag=0;

			for(i=1;i<notu;i++)
				if(abs(treer[i][j1]-treer[i][j2])!=1)
					flag=1;

			if (flag==0)
				{
				br_to_remove=j1;
				break;
				}
			}

		if (br_to_remove!=-1) break;
		}

	if (br_to_remove!=-1) break;
	}

/* writing of treenr */

for(j=0;j<br_to_remove;j++)
	{
	for(i=1;i<notu+1;i++)
		{
		treenr[i-1][j]=treer[i][j];
		}
	}

for(j=br_to_remove+1;j<notu-2;j++)
	{
	for(i=1;i<notu+1;i++)
		{
		treenr[i-1][j-1]=treer[i][j];
		}
	}

return 0;
}

/*****************************************************************************/

void Allgap(struct sequences **lineage, int *nb, int n_lineages, char *charlist)
{
int i, j, k, l;
int lg;
int *gap, ok;
char *error_message;

lg=lineage[0][0].lg;

for (i=0; i<n_lineages+1; i++)
	for (j=0; j<nb[i]; j++)
		if (lineage[i][j].lg!=lg)
			{
			error_message=Allocation(60, sizeof(char), "error_message");
			sprintf(error_message, "Sequences must be aligned: %s and %s are of different lengths.", lineage[i][j].name, lineage[0][0].name);
			Leave_program(error_message);
			}

gap=Allocation(lg, sizeof(int), "gap");

for (i=0; i<lg; i++)
	{
	gap[i]=0;

	for (j=0; j<n_lineages+1 && !gap[i]; j++)
		for (k=0; k<nb[j] && !gap[i]; k++)
			{
			ok=0;
			l=0;
			while (charlist[l] && !ok)
				{
				if (lineage[j][k].seq[i]==charlist[l]) ok=1;
				l++;
				}

			if (!ok) gap[i]=1;
			}

	if (gap[i]) for (j=0; j<n_lineages+1; j++) for (k=0; k<nb[j]; k++) lineage[j][k].seq[i]='-';
	}

gap=De_allocation(gap);
}

/*****************************************************************************/

/*returns:
gc[0]: GC in all positions
gc[1]: GC in position 1
gc[2]: GC in position 2
gc[3]: GC in position 3 */

double *GC(struct sequences s, struct options *choices)
{
int *l;
int i, j, k;
double *gc;
int n;

if (choices->seqtype) n=4;
else n=1;

l=Allocation(n, sizeof(int), "l");
gc=Allocation(n, sizeof(double), "gc");

for (i=0; i<n; i++)
	{
	gc[i]=0;
	l[i]=0;
	}

for (i=0; i<s.lg; i++)
	if (s.seq[i]=='G' || s.seq[i]=='C')
		{
		gc[0]++;
		l[0]++;

		if (choices->seqtype)
			for (j=0; j<3; j++)
				if (((i+j+1)%3)==0)
					{
					gc[3-j]++;
					l[3-j]++;
					}
		}
	else
		if (s.seq[i]=='A' || s.seq[i]=='T')
			{
			l[0]++;
			if (choices->seqtype)
				for (j=0; j<3; j++)
					if (((i+j+1)%3)==0)
						l[3-j]++;
			}

for (i=0; i<n; i++) gc[i]=Divide(gc[i], l[i]);

l=De_allocation(l);

return gc;
}

/*****************************************************************************/

struct parameters_k NC_dist(struct sequences *to_treat, int method)
{
switch (method)
	{
	case jc:
		return Poisson(to_treat, DNA_POISSON);
	case k2:
		return Kim2(to_treat);
	}

Leave_program("Error in distance method choice.");
}

/*****************************************************************************/

double NC_cov(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34, int method)
{
switch (method)
	{
	case jc:
		return CovPoisson(p12, p13, p14, p23, p24, p34, DNA_POISSON);
	case k2:
		return CovK2(p12, p13, p14, p23, p24, p34);
	}

Leave_program("Error in distance method choice.");
}

/*****************************************************************************/

struct parameters_k Kim2(struct sequences *to_treat)
{
struct parameters_k result;
double p=0., q=0.;
double a, b;
int i;

result.l=to_treat[0].lg;

for (i=0; i<to_treat[0].lg; i++)
	if (to_treat[0].seq[i]!=to_treat[1].seq[i])
		switch (to_treat[0].seq[i])
			{
			case 'A':
				if (to_treat[1].seq[i]=='G')
					{
					p++;
					break;
					}
			case 'G':
				if (to_treat[1].seq[i]=='A') p++;
				else q++;
				break;
			case 'T':
			case 'U':
				if (to_treat[1].seq[i]=='C')
					{
					p++;
					break;
					}
			case 'C':
				if (to_treat[1].seq[i]=='T' || to_treat[1].seq[i]=='U') p++;
				else q++;
				break;
			default:
				result.l--;
			}

p=Divide(p, (double)result.l);
q=Divide(q, (double)result.l);

if ((1.-2.*q)<=0 || (1.-2.*p-q)<=0) result.k=result.vk=-1;
else
	{
	result.k= -0.5 * (double)log( (1.-2.*p-q) * sqrt(1.-2.*q) );

	result.a= (double)log(1.-2.*q)/4 - (double)log(1.-2.*p-q)/2;
	result.b=-(double)log(1.-2.*q)/2;

	a=Divide(1., 1.-2.*p-q);
	b=( a + Divide(1., 1.-2.*q) )/2.;

	result.vk=Divide(a*a*p + b*b*q - (a*p+b*q)*(a*p+b*q), (double)result.l);
	}

return result;
}

/*****************************************************************************/

void *Allocation(int nb, size_t Size, char *variable_name)
{
void *result;

nb++;

if ( (result=calloc((size_t)nb, Size))!=NULL) return result;

fprintf(stderr, "Not enough memory in RRTree for %s\n", variable_name);
fprintf(stderr, "\a");
fflush(stderr);
exit(1);
}

/*****************************************************************************/

void *De_allocation(void *p)
{
if (p) free(p);
return NULL;
}

/*****************************************************************************/

int Answer(int default_answer)
{
char read[20];
int result;

Read_line(read, 20, stdin);
if (read[0]==0 || read[0]=='\n') result=default_answer;
else sscanf(read, "%d", &result);

if (result<0) Leave_program("Your choice, good bye.");

return result;
}

/*****************************************************************************/

FILE *Open_file(char mode[], char *name)
{
FILE *f;
int replace;
char *answer;

answer=Allocation(30, sizeof(char), "answer");

if (name==NULL)
	{
	name=Allocation(40, sizeof(char), "name");
	Read_line(answer, strlen(answer), stdin);
	if (answer[0]==0 || answer[0]=='\n')
		{
		answer=De_allocation(answer);
		if (mode[0]=='r' || mode[0]=='a') return stdin;
		else return stdout;
		}
	sscanf(answer, "%s", name);
	}

if (mode[0]=='w' && fopen(name, "r")!=NULL)
	{
	printf("%s exists, replace it (1/[0]) ? ", name);
	Read_line(answer, strlen(answer), stdin);

	if (answer[0]==0 || answer[0]=='\n') replace=0;
	else sscanf(answer, "%d", &replace);

	if (!replace)
		{
		do
			{
			printf("new name: ");
			Read_line(answer, strlen(answer), stdin);
			}
		while (answer[0]==0 || answer[0]=='\n');

		sscanf(answer, "%s", name);
		f=Open_file(mode, name);
		answer=De_allocation(answer);
		return f;
		}
	}

if ((f=fopen(name, mode))==NULL)
	{
	printf("file %s does not exist, input new name:\n", name);
	f=Open_file(mode, NULL);
	}

answer=De_allocation(answer);
return f;
}

/*****************************************************************************/

/*
if nb<0, read until next ";;"

BEWARE, in that case the table of sequences should have been provided (allocated)
big enough, no parachute is provided.
*/

int Read_mase(FILE *f, int nb, struct sequences *s, int choice_com)
{
char *line, *bases, *comments, *tmp;
int i, j;
int lcom=0, lline, lgbases, lgcom;
int n=0;
int the_end=0;
FILE *ff;	/* ff keeps in memory the previous position of f */

ff=Allocation(1, sizeof(FILE), "ff");
*ff=*f;

line=Allocation(MAXLEN, sizeof(char), "line");

lgbases=MAXLEN;

bases=Allocation(lgbases, sizeof(char), "bases");

lgcom=MAXLEN;

comments=Allocation(lgcom, sizeof(char), "comments");

Read_line(line, MAXLEN, f);

lline=strlen(line);

bases[0]=0;
comments[0]=0;

while (!the_end)
	{
	s[n].lg=0;
	s[n].com=NULL;
	if (n!=0 && nb<0 && line[1]==';') break;

	while (line[0]==';')
		{
		if (choice_com && line[1]!=0 && line[1]!='\n')
			{
			for (i=1; i<lline; i++) line[i-1]=line[i]; /*removes the first;*/
			lline--;
			line[lline]=0;	/*removes the \n*/
			lcom+=lline;

			if (lcom>lgcom)
				{
				lgcom+=MAXLEN;
		
				tmp=comments;
				comments=Allocation(lgcom, sizeof(char), "comments");
				memcpy(comments, tmp, lgcom-MAXLEN);
				tmp=De_allocation(tmp);

				comments[lcom]=0;
				}
			strcat(comments, line);
			}
		*ff=*f;
		Read_line(line, MAXLEN, f);
		lline=strlen(line);
		}

	if (choice_com && lcom)
		{
		if (comments[lcom-1]!='\n')
			{
			strcat(comments, "\n");
			lcom++;
			}

		s[n].com=Allocation(lcom, sizeof(char), "s[n].com");

		memcpy(s[n].com, comments, lcom);
		s[n].com[lcom]=0;

		comments[0]=0;
		lcom=0;
		}

	i=0;
	while (line[i]!='\n') i++;
	line[i]=0;			/*removes \n*/
	lline=strlen(line);

	s[n].name=Allocation(lline, sizeof(char), "s[n].name");

	memcpy(s[n].name, line, lline);

	s[n].name[lline]=0;
	
	*ff=*f;

	while (Read_line(line, MAXLEN, f)!=NULL)
		{
		if (line[0]==';')
			{
			the_end=-1;
			lline=strlen(line);
			break;
			}

		i=0;
		while (line[i]!='\n') i++;
		line[i]=0;			/*remove the \n*/
		lline=strlen(line);

		if ((s[n].lg+lline)>lgbases)
			{
			lgbases+=MAXLEN;

			tmp=Allocation(lgbases, sizeof(char), "tmp");
			memcpy(tmp, bases, lgbases-MAXLEN);
			bases=De_allocation(bases);
			bases=tmp;

			bases[s[n].lg]=0;

			strcat(bases, line);
			s[n].lg+=lline;
			}
		else
			{
			s[n].lg+=lline;
			strcat(bases, line);
			}
		
		*ff=*f;
		}

	the_end++;

	s[n].seq=Allocation(s[n].lg, sizeof(char), "s[n].seq");

	memcpy(s[n].seq, bases, s[n].lg);

	bases=De_allocation(bases);
	lgbases=MAXLEN;
	bases=Allocation(lgbases, sizeof(char), "bases");

	if ((++n)==nb) break;
	}

if (fgetc(f)!=EOF)
	*f=*ff;		/* this was suggested to me by Tal Pupko */
else	the_end=1;

line=De_allocation(line);
bases=De_allocation(bases);
comments=De_allocation(comments);

return the_end;
}

/*****************************************************************************/

struct sequences Nogap_cds(struct sequences in, char *charlist)
{
struct sequences out;
int i, j, k, l;
int ok;
char *error_message;

if ((in.lg%3)!=0)
	{
	error_message=Allocation(40+strlen(in.name), sizeof(char), "error_message");
	sprintf(error_message, "%s of length %d non divisible by 3.", in.name, in.lg);
	Leave_program(error_message);
	}

out.name=in.name;
out.com=in.com;
out.lg=in.lg;

for (i=0; i<(in.lg/3); i++)
	{
	for (j=0; j<3; j++)
		{
		ok=0;
		k=0;
		while (charlist[k])
			{
			if ((in.seq[3*i+j]==charlist[k]))
				{
				ok=1;
				break;
				}
			k++;
			}

		if (!ok)
			{
			out.lg-=3;
			break;
			}
		}
	}

out.seq=Allocation(out.lg, sizeof(char), "out.seq");

l=0;

for (i=0; i<(in.lg/3); i++)
	{
	for (j=0; j<3; j++)
		{
		ok=0;
		k=0;
		while (charlist[k])
			{
			if ((in.seq[3*i+j]==charlist[k]))
				{
				ok=1;
				break;
				}
			k++;
			}

		if (!ok) break;
		}

	if (ok)
		for (j=0; j<3; j++)
			{
			out.seq[l]=in.seq[3*i+j];
			l++;
			}
	}

return out;
}

/*****************************************************************************/

int Is_text_zero(char *text, int start)
{
int i=0;

if (!isdigit(text[start])) return 0;

while (text[start+i])
	{
	if (isdigit(text[start+i]) && text[start+i]!='0') return 0;
	if (!isdigit(text[start+i]) && text[start+i]!='.') return 1;
	i++;
	}

return 1;
}

/*****************************************************************************/

struct sequences Nogap(struct sequences in, char *charlist)
{
struct sequences out;
int i, j, k;
int ok;

out.name=in.name;
out.com=in.com;
out.lg=in.lg;

for (i=0; i<in.lg; i++)
	{
	ok=0;

	j=0;
	while (charlist[j])
		{
		if (in.seq[i]==charlist[j])
			{
			ok=1;
			break;
			}
		j++;
		}

	if (!ok) out.lg--;
	}

out.seq=Allocation(out.lg, sizeof(char), "out.seq");

k=0;

for (i=0; i<in.lg; i++)
	{
	ok=0;

	j=0;
	while (charlist[j])
		{
		if (in.seq[i]==charlist[j])
			{
			ok=1;
			break;
			}
		j++;
		}

	if (ok)
		{
		out.seq[k]=in.seq[i];
		k++;
		}
	}

return out;
}

/*****************************************************************************/

double CovPoisson(struct parameters_k p12, struct parameters_k p13, struct parameters_k p14, struct parameters_k p23, struct parameters_k p24, struct parameters_k p34, double fraction)
{
double l, k, p;	/*parameters of branch 00' of the tree*/

l=(p13.l + p14.l + p23.l + p24.l)/4.;

k=(p13.k + p14.k + p23.k + p24.k - 2*p12.k - 2*p34.k)/4.;

p=fraction*(1 - exp(-k/fraction));

if ((1-p/fraction)<=0) return -1.;
else return Divide(Divide(p*(1-p), l), (1-p/fraction)*(1-p/fraction));	/* there was a bug here identified thanks to Sandrine Hugues 16/02/2000 */
}

/*****************************************************************************/

struct parameters_k Poisson(struct sequences *to_treat, double fraction)
{
struct parameters_k result;
double p=0.;
int i;
char *error_message;

result.l=0;

if (to_treat[0].lg!=to_treat[1].lg)
	{
	error_message=Allocation(40+strlen(to_treat[0].name)+strlen(to_treat[1].name), sizeof(char), "error_message");
	sprintf(error_message, "sequences %s and %s of different length.", to_treat[0].name, to_treat[1].name);
	Leave_program(error_message);
	}

for (i=0; i<to_treat[0].lg; i++)
	if (to_treat[0].seq[i]!='-' && to_treat[1].seq[i]!='-')
		{
		result.l++;
		if(to_treat[0].seq[i]!=to_treat[1].seq[i]) p++;
		}

p=Divide(p, (double)result.l);

if ((1.-p/fraction)<=0)
	result.k=result.vk=-1.;
else
	{
	result.k=-fraction*log(1.-p/fraction);
	result.vk=Divide(Divide(p*(1.-p), result.l), (1-p/fraction)*(1-p/fraction));
	}

return result;
}

/*****************************************************************************/

int Read_fasta_gde(FILE *f, int nb, struct sequences *s, char before_name)
{
char line[MAXLEN], bases[MAXLEN];
int i, j=0;
int end=0;
FILE *ff;	/* ff keeps in memory the previous position of f */

ff=Allocation(1, sizeof(FILE), "ff");
*ff=*f;

Read_line(line, MAXLEN, f);

for (i=0; i<nb; i++)
	{
	s[i].name=Allocation(strlen(line), sizeof(char), "s[i].name");

	sscanf(line, "%*c%s", s[i].name);

	s[i].lg=0;
	
	*ff=*f;
	Read_line(line, MAXLEN, f);

	while (line[0]!=before_name)
		{
		j=0;
		while (line[j]!='\n')
			bases[s[i].lg++]=line[j++];
		*ff=*f;
		if (Read_line(line, MAXLEN, f)==NULL) break;
		}

	s[i].seq=Allocation(s[i].lg, sizeof(char), "s[i].seq");

	for (j=0; j<s[i].lg; j++)
		s[i].seq[j]=bases[j];
	s[i].seq[s[i].lg]=0;

	s[i].com=NULL;
	}

if (fgetc(f)!=EOF)
	*f=*ff;		/* this was suggested to me by Tal Pupko */
else	end=1;

return end;
}

/*****************************************************************************/

int Read_lintre(FILE *f, int nb, struct sequences *s)
{
char one_char, name[NAMELEN];
int i, j, end=0;
long file_position;
char pc_char='\r';


for (i=0; i<nb; i++)
	{
	fscanf(f, "%s", name);
	s[i].name=Allocation(strlen(name)+1, sizeof(char), "s[i].name");
	strcpy(s[i].name, name);
	s[i].name[strlen(name)]=0;

	file_position=ftell(f);	/* mark begining of sequence */

	s[i].lg=0;
	one_char=' ';
	while (one_char && one_char!='\n' && one_char!=EOF && one_char!=pc_char)	/* count length of sequence */
		{
		one_char=fgetc(f);
		s[i].lg++;
		}

	s[i].seq=Allocation(s[i].lg+1, sizeof(char), "s[i].seq");

	fseek(f, file_position, 0);	/* return to begining of sequence */

	Read_line(s[i].seq, s[i].lg, f);

	one_char=fgetc(f);
	while (one_char && one_char!='\n' && one_char!=EOF && one_char!=pc_char) one_char=fgetc(f);	/* go beyond '\n' */

	for (j=0; j<s[i].lg; j++)
		if (s[i].seq[j]=='?') s[i].seq[j]='-';
	s[i].seq[s[i].lg]=0;

	s[i].com=NULL;

	file_position=ftell(f);	/* mark end of reading */
	one_char=fgetc(f);
	while (one_char=='\n' || one_char=='\t' || one_char==' ' || one_char==pc_char)	one_char=fgetc(f);
	if (one_char==EOF) end=1;
	fseek(f, file_position, 0);	/* return to end of reading */
	}

return end;
}

/*****************************************************************************/

int Read_phylip(FILE *f, int to_read, int read, struct sequences *s, int lg_name)	/* reads the interleaved PHYLIP format */
{
char line[MAXLEN];
int lg, nb;
int i, j;
char *error_message;

Get_full_line(line, MAXLEN, f);
sscanf(line, "%d%d", &nb, &lg);

nb = nb - read - to_read;

if (nb<0)
	{
	error_message=Allocation(70, sizeof(char), "error_message");
	sprintf(error_message, "impossible to read %d+%d sequences, file only contains %d.", read, to_read, nb+read+to_read);
	Leave_program(error_message);
	}

for (i=0; i<read; i++) Get_full_line(line, MAXLEN, f);

for (i=0; i<to_read; i++)
	{
	Get_full_line(line, MAXLEN, f);

	if (lg_name)
		{
		s[i].name=Allocation(lg_name+1, sizeof(char), "s[i].name");

		for (j=0; j<lg_name; j++)
			{
			if (line[j]==' ') break;
			s[i].name[j]=line[j];
			}
		s[i].name[j]=0;
		}
	else
		{
		while (line[lg_name] && line[lg_name]!=' ' && line[lg_name]!='\t' && line[lg_name]!='\n') lg_name++;

		s[i].name=Allocation(lg_name, sizeof(char), "s[i].name");

		for (j=0; j<lg_name; j++) s[i].name[j]=line[j];
		s[i].name[j]=0;

		lg_name=0;
		}

	s[i].lg=0;
	s[i].com=NULL;
	s[i].seq=Allocation(lg+1, sizeof(char), "s[i].seq");

	if (lg_name) j=lg_name;
	else j=strlen(s[i].name);

	while (line[j] && line[j]!='\n')
		{
		if (line[j]!=' ' && line[j]!='\t')
			{
			s[i].seq[s[i].lg]=line[j];
			s[i].lg++;
			}
		j++;
		}
	}

for (i=0; i<nb; i++) Get_full_line(line, MAXLEN, f);

while (s[0].lg<lg)
	{
	for (i=0; i<read; i++) Get_full_line(line, MAXLEN, f);

	for (i=0; i<to_read; i++)
		{
		Get_full_line(line, MAXLEN, f);

		j=0;

		while (line[j] && line[j]==' ' || line[j]=='\t') j++;

		while (line[j] && line[j]!='\n')
			{
			if (line[j]!=' ' && line[j]!='\t')
				{
				s[i].seq[s[i].lg]=line[j];
				s[i].lg++;
				}
			j++;
			}
		}

	for (i=0; i<nb; i++) Get_full_line(line, MAXLEN, f);
	}

for (i=0; i<to_read; i++) s[i].seq[s[i].lg]=0;

rewind(f);

if (nb) return 0;
else return 1;	/* end of file */
}

/*****************************************************************************/

int Read_mega(FILE *f, int to_read, int read, struct sequences *s)
{
char line[MAXLEN];
int i, j;
char *error_message, *tmp, *go_on;
int not_read=0;

go_on=Read_line(line, MAXLEN, f);	/* "#mega" line */
go_on=Read_line(line, MAXLEN, f);	/* title line */
while (line[0]!='#') go_on=Read_line(line, MAXLEN, f);	/* go to the first sequence */

for (i=0; i<read; i++) go_on=Read_line(line, MAXLEN, f);	/* jump the sequences previously read */

for (i=0; i<to_read && go_on; i++)
	{
	j=0;
	while (line[j] && line[j]!=' ' && line[j]!='\t' && line[j]!='\n') j++;
	s[i].name=Allocation(j, sizeof(char), "s[i].name");
	j=0;
	while (line[j+1] && line[j+1]!=' ' && line[j+1]!='\t')
		{
		s[i].name[j]=line[j+1];	/* +1 because of the initial '#' */
		j++;
		}
	s[i].name[j]=0;
	while (line[j]==' ' || line[j]=='\t') j++;

	s[i].lg=0;
	s[i].com=NULL;
	s[i].seq=Allocation(strlen(line)-j, sizeof(char), "s[i].seq");

	while (line[j] && line[j]!='\n')
		switch (line[j])
			{
			case ' ':
			case '\t':
				j++;
				break;
			case '"':
				j++;
				while (line[j] && line[j]!='\n' && line[j]!='"') j++;
				j++; 	/* jump the '"' itself */
				break;
			default:
				s[i].seq[s[i].lg]=line[j];
				s[i].lg++;
				j++;
			}

	go_on=Read_line(line, MAXLEN, f);
	}

while (go_on && line[0]!='\n')
	{
	go_on=Read_line(line, MAXLEN, f);	/* jump the sequences we're not reading now */
	not_read++;
	}

while (go_on && line[0]=='\n')
	go_on=Read_line(line, MAXLEN, f);	/* jump the sequences we're not reading now */

while (go_on)
	{
	for (i=0; i<read; i++) go_on=Read_line(line, MAXLEN, f);	/* jump the sequences previously read */

	for (i=0; i<to_read; i++)
		{
		j=0;
		while (line[j]!=' ' && line[j]!='\t') j++;	/* skip name */
		while (line[j]==' ' || line[j]=='\t') j++;

		tmp=Allocation(s[i].lg+1, sizeof(char), "tmp");
		strcpy(tmp, s[i].seq);
		s[i].seq=De_allocation(s[i].seq);
		s[i].seq=Allocation(s[i].lg+strlen(line)-j, sizeof(char), "s[i].seq");
		strcpy(s[i].seq, tmp);
		tmp=De_allocation(tmp);

		while (line[j] && line[j]!='\n')
			switch (line[j])
				{
				case ' ':
				case '\t':
					j++;
					break;
				case '"':
					j++;
					while (line[j] && line[j]!='\n' && line[j]!='"') j++;
					break;
				default:
					s[i].seq[s[i].lg]=line[j];
					s[i].lg++;
					j++;
				}

		go_on=Read_line(line, MAXLEN, f);
		}

	while (go_on && line[0]!='\n') go_on=Read_line(line, MAXLEN, f);	/* jump the sequences we're not reading now */
	while (go_on && line[0]=='\n') go_on=Read_line(line, MAXLEN, f);
	}

for (i=0; i<to_read; i++) s[i].seq[s[i].lg]=0;

rewind(f);

return 1-not_read;	/* end of file? */
}

/*****************************************************************************/

int Read_phyltest(FILE *f, int to_read, int read, struct sequences *s)
{
int i, j, k;
int nb, lg;
int go_on=1;
char *line, *error_message;

line=Allocation(MAXLEN, sizeof(char), "line");

Get_full_line(line, MAXLEN, f);	/* description line */
Get_full_line(line, MAXLEN, f);	/* data type */
Get_full_line(line, MAXLEN, f);	/* sequence number and length */

sscanf(line, "%d%d", &nb, &lg);
nb = nb - read - to_read;
if (nb<0)
	{
	error_message=Allocation(70, sizeof(char), "error_message");
	sprintf(error_message, "impossible to read %d+%d sequences, file only contains %d.", read, to_read, nb+read+to_read);
	Leave_program(error_message);
	}

for (i=0; i<to_read; i++)
	{
	s[i].lg=lg;
	s[i].seq=Allocation(lg+1, sizeof(char), "s[i].seq");
	s[i].com=NULL;
	}

if (lg>MAXLEN)
	{
	line=De_allocation(line);
	line=Allocation(lg+1, sizeof(char), "line");
	}
else	lg=MAXLEN;	/* maximum length of lines for this file */

nb=0;

while (line[0]!='#' && go_on) go_on=Get_full_line(line, lg, f);	/* going to the first sequence */

while (nb<read)
	{
	while (line[0]=='#' && go_on) go_on=Get_full_line(line, lg, f);	/* skipping sequences previously read */
	while (line[0]!='#' && go_on) go_on=Get_full_line(line, lg, f);
	nb++;
	}

for (i=0; i<to_read; i++)
	{
	while (line[0]!='#' && go_on) go_on=Get_full_line(line, lg, f);
	j=0;
	while (line[j+1] && line[j+1]!='{') j++;	/* '{' starts the lineage name, so ends the sequence name */
	s[i].name=Allocation(j+1, sizeof(char), "s[i].name");
	j=0;
	while (line[j+1] && line[j+1]!='{')
		{
		s[i].name[j]=line[j+1];	/* +1 because of the first '#' */
		j++;
		}
	if (s[i].name[j-1]=='_') s[i].name[j-1]=0;
	else s[i].name[j]=0;

	k=0;
	go_on=Get_full_line(line, lg, f);
	while (line[0]!='#' && go_on)
		{
		j=0;
		while (line[j])
			{
			if (line[j]!=' ' && line[j]!='\t' && line[j]!='\n')
				{
				s[i].seq[k]=line[j];
				k++;
				}
			j++;
			}

		go_on=Get_full_line(line, lg, f);
		}
	}

rewind(f);

return 1-go_on;	/* 1 gives 0, 0 gives 1 */
}

/*****************************************************************************/

int Get_full_line(char *line, int maxlen, FILE *f)
{
int i;

while (Read_line(line, maxlen, f))
	while (line[0])
		{
		if (line[0]!=' ' && line[0]!='\t' && line[0]!='\n') return 1;

		i=1;
		while (line[i])
			{
			line[i-1]=line[i];
			i++;
			}
		line[i-1]=0;
		}

return 0;
}

/*****************************************************************************/

int What_file_format(FILE *f, FILE *in_command_file)
{
char first_char=' ', *read, line[MAXLEN];
int result, i;

while (first_char==' ' || first_char=='\t' || first_char=='\n') first_char=fgetc(f);

rewind(f);

read=Read_in_command_file_string(in_command_file, "format");
if (read)
	{
	i=0;
	while (read[i])
		{
		read[i]=toupper(read[i]);
		i++;
		}
	}

switch (first_char)
	{
	case ';':
		if (read && strcmp(read, "MASE")==0) return mase;

		printf("MASE format (0/[1]) ? ");
		result=Answer(1);
		if (result)  return mase;
		break;
	case '>':
		if (read && strcmp(read, "FASTA")==0) return fasta;

		printf("FASTA format (0/[1]) ? ");
		result=Answer(1);
		if (result)  return fasta;
		break;
	case 'C':
		if (read && strcmp(read, "CLUSTAL")==0) return clustal;

		printf("CLUSTAL format (0/[1]) ? ");
		result=Answer(1);
		if (result)  return clustal;
		break;
	case '%':
		if (read && strcmp(read, "GDE")==0) return gde;

		printf("GDE format (0/[1]) ? ");
		result=Answer(1);
		if (result)  return gde;
		break;
	case '#':
		Read_line(line, MAXLEN, f);
		rewind(f);
		i=0;
		while (line[i])
			{
			line[i]=toupper(line[i]);
			i++;
			}
		if (strstr(line, "MEGA"))
			{
			if (read && strcmp(read, "MEGA")==0) return mega;

			printf("interleaved MEGA format (0/[1]) ? ");
			result=Answer(1);
			if (result)  return mega;
			}
		else
			{
			if (read && strcmp(read, "NEXUS")==0) return nexus;

			printf("NEXUS format (0/[1]) ? ");
			result=Answer(1);
			if (result)  return nexus;
			}
		break;
	default:	/* these are not recognisable by their first character */
		for (i=0; i<4; i++) Get_full_line(line, MAXLEN, f);	/* the 1st three lines are various text */
		rewind(f);
		if (line[0]=='#')	/* then the 1st sequence should start with a '#' */
			{
			if (read && strcmp(read, "PHYLTEST")==0) return phyltest;

			printf("PHYLTEST format (0/[1]) ? ");
			result=Answer(1);
			if (result)  return phyltest;
			}
		else
			{
			if (isdigit(first_char))
				{
				if (read && strcmp(read, "PHYLIP")==0) return phylip;

				printf("interleaved PHYLIP format (0/[1]) ? ");
				result=Answer(1);
				if (result)  return phylip;
				}
			else	/* I didn't find anything allowing to recognize LINTRE files, so it's the default. */
				{
				if (read && strcmp(read, "LINTRE")==0) return lintre;

				printf("LINTRE format (0/[1]) ? ");
				result=Answer(1);
				if (result)  return lintre;
				}
			}
	}

printf("Please specify a file format: MASE (%d), FASTA (%d), interleaved PHYLIP (%d), CLUSTAL (%d), GDE (%d), NEXUS (%d), LINTRE (%d), PHYLTEST (%d):\n", mase, fasta, phylip, clustal, gde, nexus, lintre, phyltest);
result=Answer(-1);

if (result==mase || result==fasta || result==phylip || result==clustal || result==gde || result==nexus || result==lintre || result==phyltest)
	{
	printf("Your choice of format will be used, but beware that since it was not recognized by the program, there may be problems in the lecture of the sequence data.\n");
	return result;
	}

Leave_program("Unknown file format, sorry.");
}

/*****************************************************************************/

void What_sequence_type(char *seq, struct options *choices, FILE *in_command_file)
{
int i, lg=0, stop;
double numACGT=0;
char *read;

choices->seqtype=NOT_IN_FILE;

read=Read_in_command_file_string(in_command_file, "type");

if (read)
	{
	i=0;
	while (read[i])
		{
		read[i]=toupper(read[i]);
		i++;
		}

	if (strstr(read, "NC"))
		{
		choices->seqtype=NC;
		return;
		}

	if (strstr(read, "PROT"))
		{
		choices->seqtype=PROT;
		return;
		}

	if (strstr(read, "CDS"))
		{
		choices->seqtype=CDS;
		}
	}

if (choices->seqtype==NOT_IN_FILE)
	{
	i=0;
	while (seq[i])
		{
		if (isalpha(seq[i]))
			{
			lg++;
			if (seq[i]=='A' || seq[i]=='C' || seq[i]=='G' || seq[i]=='T' || seq[i]=='U') numACGT++;
			}
		i++;
		}

	if (Divide(numACGT, lg)>=0.7)	/* if it's mostly nucleotides, it's probably DNA */
		{
		stop=lg%3;	/* if it's not dividable by 3, it's probably not a coding sequence */
	
		if (stop==0)	/* if it has stop codons in the middle, it's probably not a coding sequence */
			{
			i=0;
			while (seq[3*i+3] && stop==0)	/* don't look at the last codon */
				{
				if ((seq[3*i]=='T' || seq[3*i]=='U') && ((seq[3*i+1]=='A' && seq[3*i+2]=='G') || (seq[3*i+1]=='G' && seq[3*i+2]=='A')))
					stop++;
				i+=3;
				}
			}
		
		if (stop==0)
			{
			printf("non coding DNA sequences (0), coding DNA sequences [1], or amino-acid sequences (2) ? ");
			choices->seqtype=Answer(CDS);
			}
		else
			{
			printf("non coding DNA sequences [0], coding DNA sequences (1), or amino-acid sequences (2) ? ");
			choices->seqtype=Answer(NC);
			}
		}
	else
		{
		printf("non coding DNA sequences (0), coding DNA sequences (1), or amino-acid sequences [2] ? ");
		choices->seqtype=Answer(PROT);
		}
	}

if (choices->seqtype==CDS)
	{
	choices->code_mt=NOT_IN_FILE;

	read=Read_in_command_file_string(in_command_file, "code");
	if (read)
		{
		if (read[0]=='n' || read[0]=='N') choices->code_mt=0;	/* nuclear */
		if (read[0]=='m' || read[0]=='M') choices->code_mt=1;	/* mitochondrial */
		}

	if (choices->code_mt==NOT_IN_FILE)
		{
		printf("code: nuclear [0] or mitochondrial (1) ? ");
		choices->code_mt=Answer(0);
		}
	}
}

/*****************************************************************************/

int Read_clustal(FILE *f, int to_read, int read, struct sequences *s)
{
char line[MAXLEN], *p;
int end=0, *lg_alloc;
int i, j;

if (Read_line(line, MAXLEN, f)==NULL) end=1;
while (line[0]!='\n') if (Read_line(line, MAXLEN, f)==NULL) end=1;
while (line[0]=='\n') if (Read_line(line, MAXLEN, f)==NULL) end=1;

lg_alloc=Allocation(to_read, sizeof(int), "lg_alloc");

for (i=0; i<to_read; i++)
	{
	lg_alloc[i]=MAXLEN;

	s[i].name=Allocation(NAMELEN_CLUSTAL, sizeof(char), "s[i].name");
	s[i].name[0]=0;
	s[i].com=NULL;
	s[i].seq=Allocation(lg_alloc[i], sizeof(char), "s[i].seq");
	s[i].lg=0;
	}

while (!end)
	{
	for (i=0; i<read && !end; i++)
		if (Read_line(line, MAXLEN, f)==NULL) end=1;

	for (i=0; i<to_read && !end; i++)
		{
		if (!s[i].name[0])
			for (j=0; j<NAMELEN_CLUSTAL; j++)
				{
				if (line[j]==' ') s[i].name[j]=0;
				else s[i].name[j]=line[j];
				}
		else j=NAMELEN_CLUSTAL;

		while (line[j] && line[j]!='\n')
			{
			if (line[j]!=' ' && line[j]!='\t')
				{
				s[i].seq[s[i].lg]=line[j];
				s[i].lg++;
				
				if (s[i].lg>=(lg_alloc[i]-1))
					{
					lg_alloc[i]+=MAXLEN;

					p=s[i].seq;
					s[i].seq=Allocation(lg_alloc[i], sizeof(char), "s[i].seq");
					strcpy(s[i].seq, p);
					p=De_allocation(p);
					}
				}
				
			j++;
			}

		if (Read_line(line, MAXLEN, f)==NULL) end=1;
		}

	i=0;
	while (line[0]!='\n' && !end)
		{
		if (Read_line(line, MAXLEN, f)==NULL) end=1;
		i++;
		}
	while (line[0]=='\n' && !end) if (Read_line(line, MAXLEN, f)==NULL) end=1;
	}

rewind(f);

if (i>1) return 0;
else return 1;
}

/*****************************************************************************/

int Read_nexus(FILE *f, int to_read, int read, struct sequences *s)
{
char line[MAXLEN], *p;
int end=0, *lg_alloc, return_end;
int i, j;

/* reads one sequence too many? */

lg_alloc=Allocation(to_read, sizeof(int), "lg_alloc");

for (i=0; i<to_read; i++)
	{
	lg_alloc[i]=MAXLEN;

	s[i].lg=0;
	s[i].com=NULL;
	s[i].name=NULL;
	s[i].seq=Allocation(lg_alloc[i], sizeof(char), "s[i].seq");
	}

line[0]=0;

while (!strstr(line, "MATRIX"))
	{
	Get_full_line(line, MAXLEN, f);
	/* put in capitals */
	i=0;
	while (line[i])
		{
		line[i]=toupper(line[i]);
		i++;
		}
	}

i=0;

if (Get_full_line(line, MAXLEN, f)==0) end=1;

while (!end && line[0]!=';')
	{
	if (i>=to_read || line[0]=='\n')
		{
		j=0;

		while (line[0]!='\n' && line[0]!=';')
			{
			if (Get_full_line(line, MAXLEN, f)==0) end=1;
			j++;
			}

		while (line[0]=='\n' && line[0]!=';')
			{
			if (Get_full_line(line, MAXLEN, f)==0) end=1;
			j++;
			}

		i=0;

		if (j<=1 && line[0]==';') return_end=1;
		else return_end=0;
		}
	else
		{
		if (i==0)
			for (j=0; j<read; j++)
				if (Get_full_line(line, MAXLEN, f)==0 || line[0]==';') end=1;

		if (!s[i].name)
			{
			j=0;
			while (line[j]!=' ' && line[j]!='\t') j++;

			s[i].name=Allocation(j+1, sizeof(char), "s[i].name");
		
			j=0;
		
			while (line[j]!=' ' && line[j]!='\t')
				{
				s[i].name[j]=line[j];
				j++;
				}
			
			s[i].name[j]=0;
			}
		else
			{
			j=0;
			while (line[j]!=' ' && line[j]!='\t') j++;
			}

		while (line[j]==' ' || line[j]=='\t') j++;

		if ((s[i].lg+strlen(line)-j)>=lg_alloc[i])
			{
			lg_alloc[i]+=MAXLEN;
			p=s[i].seq;
			s[i].seq=Allocation(lg_alloc[i], sizeof(char), "s[i].seq");
			strcpy(s[i].seq, p);
			p=De_allocation(p);
			}

		while (line[j] && line[j]!='\n')
			{
			if (line[j]!=' ' && line[j]!='\t')
				{
				s[i].seq[s[i].lg]=line[j];
				s[i].lg++;
				}

			j++;
			}

		i++;

		if (Get_full_line(line, MAXLEN, f)==0 || line[0]==';') end=1;
		}
	}

rewind(f);

return return_end;
}


/*****************************************************************************/

char *Unsignificant_branches_to_zero(char *tree_in, int limit)
{
char *tree_out, boot[10], *p;
int i, j, k, l, i_out;
int lg, val_boot;
char car;
int already_a_distance=0;

lg=strlen(tree_in)+1;

tree_out=Allocation(lg, sizeof(char), "tree_out");

i=0;
i_out=0;
while (tree_in[i])
	{
	if (tree_in[i-1]==')' && isdigit(tree_in[i]))
		{
		j=i;

		while (isdigit(tree_in[j]))
			{
			boot[j-i]=tree_in[j];
			j++;
			}
		boot[j-i]=0;

		val_boot=atoi(boot);

		lg+=8;
		p=tree_out;
		tree_out=Allocation(lg, sizeof(char), "tree_out");
		strcpy(tree_out, p);
		p=De_allocation(p);

		tree_out[i_out]=tree_in[i];
		i_out++;

		i=j;
		while (tree_in[i]!=')' && tree_in[i]!=',') i++;

		tree_out[i_out]=':';
		i_out++;
		tree_out[i_out]='0';
		i_out++;
		tree_out[i_out]='.';
		i_out++;

		if (val_boot<limit)
			{
			for (k=0; k<4; k++) tree_out[i_out+k]='0';
			i_out+=k;
			}
		else
			{
			for (k=0; k<4; k++) tree_out[i_out+k]='1';
			i_out+=k;
			}

		already_a_distance=1;
		}
	else
		{
		if (i && (tree_in[i]==',' || tree_in[i]==')') && !already_a_distance)
			{
			j=i;
			while (j && tree_in[j]!=':') j--;

			if (tree_in[j]!=':')
				{
				lg+=8;
				p=tree_out;
				tree_out=Allocation(lg, sizeof(char), "tree_out");
				strcpy(tree_out, p);
				p=De_allocation(p);

				tree_out[i_out]=':';
				i_out++;
				tree_out[i_out]='0';
				i_out++;
				tree_out[i_out]='.';
				i_out++;
				for (k=0; k<4; k++) tree_out[i_out+k]='1';
				i_out+=k;
				}

			tree_out[i_out]=tree_in[i];
			i_out++;
			i++;

			already_a_distance=1;
			}
		else
			{
			tree_out[i_out]=tree_in[i];
			i_out++;
			i++;

			already_a_distance=0;
			}
		}
	}

tree_in=De_allocation(tree_in);

return tree_out;
}

/*****************************************************************************/

char *Make_unsignificant_dichotomies(char *tree_in, int previous_nb_add)
{
int i, j, k, l, m;
int nb_coma=0, nb_parenthesis=0, nb_add, lg_out;
char *tree_out;
int add_parenthesis;

i=0;
while (tree_in[i]!=';')
	{
	switch (tree_in[i])
		{
		case ',':
			nb_coma++;
			break;
		case '(':
			nb_parenthesis++;
			break;
		}

	i++;
	if (!tree_in[i]) Leave_program("error in tree format, must end with ';'");
	}

nb_add=nb_coma - nb_parenthesis;

if (nb_add==0 || nb_add==1) return tree_in;	/* 0 for a rooted tree, 1 for an unrooted tree */

if (nb_add<0) Leave_program("error in tree format, more '(' than ','");

lg_out=strlen(tree_in) + nb_add*9;
tree_out=Allocation(lg_out, sizeof(char), "tree_out");

i=0;
while (tree_in[i])
	{
	tree_out[i]=tree_in[i];
	i++;
	}
while (i<lg_out)
	{
	tree_out[i]=0;
	i++;
	}

i=1;
while (nb_add && tree_out[i]!=';')
	{
	add_parenthesis=0;

	if (tree_out[i]=='(')
		{
		j=1;
		nb_parenthesis=1;
		nb_coma=0;

		while (nb_parenthesis && tree_out[i+j])
			{
			switch (tree_out[i+j])
				{
				case ')':
					nb_parenthesis--;
					nb_coma--;
					break;
				case '(':
					nb_parenthesis++;
					break;
				case ',':
					nb_coma++;
				}

			j++;
			}

		if (nb_coma>0)
			add_parenthesis=1;
		else
			{
			while (tree_out[i+j] && tree_out[i+j]!='(')
				{
				if (tree_out[i+j]==',') nb_coma++;
				j++;
				}

			if (nb_coma>1) add_parenthesis=1;
			}
		}

	if (add_parenthesis)
		{
		for (j=strlen(tree_out); j>i; j--) tree_out[j]=tree_out[j-1];

		j=1;
		add_parenthesis=0;

		while (!add_parenthesis)
			{
			while (tree_out[i+j]!=';' && tree_out[i+j]!=',') j++;
			if (tree_out[i+j]!=';') j++;
			while (tree_out[i+j]!=';' && tree_out[i+j]!=',' && tree_out[i+j]!=')' && tree_out[i+j]!='(') j++;
			if (tree_out[i+j]==',' || tree_out[i+j]==';') add_parenthesis=1;
			}

		for (k=strlen(tree_out)-1; k>=i+j; k--) tree_out[k+8]=tree_out[k];

		tree_out[i+j]=':';
		j++;
		tree_out[i+j]='0';
		j++;
		tree_out[i+j]='.';
		j++;
		tree_out[i+j]='0';
		j++;
		tree_out[i+j]='0';
		j++;
		tree_out[i+j]='0';
		j++;
		tree_out[i+j]='0';
		j++;
		tree_out[i+j]=')';

		nb_add--;
		i++;
		}

	i++;
	}

nb_coma=0;
nb_parenthesis=0;
i=0;
while (tree_out[i]!=';')
	{
	switch (tree_out[i])
		{
		case ',':
			nb_coma++;
			break;
		case '(':
			nb_parenthesis++;
			break;
		}

	i++;
	}

if (nb_coma - nb_parenthesis == previous_nb_add) return tree_in;

if (nb_coma>nb_parenthesis) return Make_unsignificant_dichotomies(tree_out, nb_add);

tree_in=De_allocation(tree_in);

return tree_out;
}

/*****************************************************************************/

char *Read_tree(FILE *ftree, struct options *choices)
{
int i;
char one_char, *tree_parenthesis;

rewind(ftree);
i=0;

one_char=fgetc(ftree);
if (one_char=='[')
	while (one_char!=']')
		one_char=fgetc(ftree);
while (one_char!=';')
	{
	if ((one_char=fgetc(ftree))==EOF) Leave_program("uncorrect tree file format: ';' missing");
	i++;
	}
tree_parenthesis=Allocation(i+2, sizeof(char), "tree_parenthesis");

rewind(ftree);
i=0;

one_char=fgetc(ftree);

while (one_char!='(')
	{
	if (one_char=='[')
		while (one_char!=']') one_char=fgetc(ftree);
	one_char=fgetc(ftree);
	}

while (one_char!=';')
	{
	if (one_char!='\n' && one_char!=' ' && one_char!='\t')
		{
		tree_parenthesis[i]=one_char;
		i++;
		}
	one_char=fgetc(ftree);
	}
tree_parenthesis[i]=one_char;
i++;
tree_parenthesis[i]=0;

tree_parenthesis=Make_unsignificant_dichotomies(tree_parenthesis, -2);

if (choices->support_limit) tree_parenthesis=Unsignificant_branches_to_zero(tree_parenthesis, choices->support_limit);

return tree_parenthesis;
}

/*****************************************************************************/

void Leave_program(char *message)
{
fprintf(stderr, "\a");	/* bip! */
fprintf(stderr, "%s\n", message);
fflush(stderr);

printf("Premature end of program\n");
fflush(stdout);

exit(1);
}

/*****************************************************************************/

double Exact_probability(double obs_value)
{
double Integral(double a, double b);
double old_result, new_result=0;
double a, b;
int nb=1;

obs_value=fabs(obs_value);

old_result=1.+PRECISION;

while (fabs(old_result-new_result)>=PRECISION)
	{
	old_result=new_result;
	new_result=0.;

	a=0.;
	b=Divide(obs_value, (double)nb);

	while (b<=obs_value)
		{
		new_result+=Integral(a, b);
		a+=Divide(obs_value, nb);
		b+=Divide(obs_value, nb);
		}

	nb++;
	}

new_result=1.-2*new_result;

if (new_result>PRECISION) return new_result;
else return PRECISION;
}

/********************************************/

double Integral(double a, double b)
{
double fa, fb;

fa=1./sqrt(2.*PI_CST)*exp(-0.5*a*a);
fb=1./sqrt(2.*PI_CST)*exp(-0.5*b*b);

return (b-a)*(fa+fb)/2.;
}

/********************************************/

/*
Reads lines of characters from a file through fgets, getting rid of
the '^M' which may occur in a file saved on a PC.
char *line must be allocated previously, at maxlen.
*/

char *Read_line(char *line, int maxlen, FILE *f)
{
char pc_char='\r';
int i, j, lg=0;

if (f==stdin)
	{
	gets(line);
	return line;
	}

if (fgets(line, maxlen, f)==NULL) return NULL;

i=0;
while (line[i])
	{
	if (line[i]==pc_char)
		{
		if (!lg) while(line[lg]) lg++;

		for (j=i; j<lg; j++)
			line[j]=line[j+1];

		lg--;
		line[lg]=0;
		}

	i++;
	}

return line;
}

/********************************************/

void Free_sequences(struct sequences *s, int nb)
{
int i;

for (i=0; i<nb; i++)
	{
	s[i].seq=De_allocation(s[i].seq);
	s[i].name=De_allocation(s[i].name);
	s[i].com=De_allocation(s[i].com);
	s[i].lg=0;
	}
}

/********************************************/

struct sequences Copy_sequence(struct sequences original)
{
struct sequences copy;
int i;

copy.lg=original.lg;

copy.seq=Allocation(original.lg+1, sizeof(char), "copy.seq");
for (i=0; i<original.lg; i++) copy.seq[i]=original.seq[i];
copy.seq[i]=0;

copy.name=Allocation(strlen(original.name)+1, sizeof(char), "copy.name");
i=0;
while (original.name[i])
	{
	copy.name[i]=original.name[i];
	i++;
	}
copy.name[i]=0;

if (original.com)
	{
	copy.com=Allocation(strlen(original.com)+1, sizeof(char), "copy.com");
	i=0;
	while (original.com[i])
		{
		copy.com[i]=original.com[i];
		i++;
		}
	copy.com[i]=0;
	}
else
	copy.com=NULL;

return copy;
}

/********************************************/

double Divide(double a, double b)
{
if (b==0) return -1;
else return a/b;
}

/********************************************/

int Read_in_command_file_number(FILE *f, char *search)
{
char *line, *search_upper;
int i, result;

if (f==NULL) return NOT_IN_FILE;

search_upper=Allocation(strlen(search)+1, sizeof(char), "search_upper");

i=0;
while (search[i])
	{
	search_upper[i]=toupper(search[i]);
	i++;
	}

line=Allocation(MAXLEN, sizeof(char), "line");

rewind(f);

while (Read_line(line, MAXLEN, f))
	{
	i=0;
	while (line[i] && line[i]!=':')
		{
		line[i]=toupper(line[i]);
		i++;
		}

	if (line[0]!='#' && strstr(line, search_upper))
		{
		i=0;
		while (line[i]!=':') i++;

		sscanf(line+i+1, "%d", &result);

		return result;
		}
	}

line=De_allocation(line);

return NOT_IN_FILE;
}

/********************************************/

char *Read_in_command_file_string(FILE *f, char *search)
{
char *line, *result, *search_upper;
int i, j;

if (f==NULL) return NULL;

search_upper=Allocation(strlen(search)+1, sizeof(char), "search_upper");

i=0;
while (search[i])
	{
	search_upper[i]=toupper(search[i]);
	i++;
	}

line=Allocation(MAXLEN, sizeof(char), "line");

rewind(f);

while (Read_line(line, MAXLEN, f))
	if (line[0]!='#' && line[0]!='\n')
		{
		i=0;
		while (line[i] && line[i]!=':')
			{
			line[i]=toupper(line[i]);
			i++;
			}

		if (strstr(line, search_upper))
			{
			i=0;
			while (line[i]!=':') i++;

			result=Allocation(strlen(line)-i, sizeof(char), "result");

			i++;	/* skip the ':' */
			while (line[i]==' ' || line[i]=='\t') i++;

			j=0;
			while (line[i] && line[i]!='\n')
				{
				result[j]=line[i];
				i++;
				j++;
				}

			line=De_allocation(line);

			return result;
			}
		}

line=De_allocation(line);

return NULL;
}

/*****************************************************************************/

int Nbseq_mase(FILE *f)
{
int nb=0;
char line[MAXLEN];
int a=1;

while (fgets(line, MAXLEN, f)!=NULL) {
	if (a && line[0]==';') {
		nb++;
		a=0;
		}
	if (!a && line[0]!=';')
		a=1;
	}

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_fasta_gde_phyltest(FILE *f, char before_name)
{
int nb=0;
char line[MAXLEN];
int a=1;

while (fgets(line, MAXLEN, f)!=NULL)
	if (line[0]==before_name) nb++;

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_phylip(FILE *f)
{
int nb;

fscanf(f, "%d", &nb);

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_clustal(FILE *f)
{
int nb=0;
char line[MAXLEN];

if (Read_line(line, MAXLEN, f)==NULL) return 0;
while (line[0]!='\n') if (Read_line(line, MAXLEN, f)==NULL) return 0;
while (line[0]=='\n') if (Read_line(line, MAXLEN, f)==NULL) return 0;

while (line[0]!='\n')
	{
	nb++;
	if (Read_line(line, MAXLEN, f)==NULL)
		{
		rewind(f);
		return nb;
		}
	}

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_nexus(FILE *f)
{
int nb=0, i;
char line[MAXLEN];

line[0]=0;
while (!strstr(line, "MATRIX"))
	{
	Get_full_line(line, MAXLEN, f);
	/* put in capitals */
	i=0;
	while (line[i])
		{
		line[i]=toupper(line[i]);
		i++;
		}
	}
Get_full_line(line, MAXLEN, f);
while (line[0]!='\n' && line[0]!=';')
	{
	nb++;
	if (Read_line(line, MAXLEN, f)==0)
		{
		rewind(f);
		return nb;
		}
	}

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_lintre(FILE *f)
{
int nb=0;
char line[MAXLEN];

while (Read_line(line, MAXLEN, f)) nb++;

rewind(f);

return nb;
}

/*****************************************************************************/

int Nbseq_mega(FILE *f)
{
int nb=0;
char line[MAXLEN];

Read_line(line, MAXLEN, f);
Read_line(line, MAXLEN, f);

while (line[0]!='#') Read_line(line, MAXLEN, f);
while (line[0]=='#')
	{
	nb++;
	if (Read_line(line, MAXLEN, f)==NULL)
		{
		rewind(f);
		return nb;
		}
	}

rewind(f);

return nb;
}

/*****************************************************************************/
/*****************************************************************************/
