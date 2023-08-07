#include <cstring>
#include <iostream>
#include <fstream>
#include "hictk/file.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

int flipSign(const char *genome, double *x, int n, char *chr, int binsize);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o observed][-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <outfile>: eigenvectors output file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
}

int main(int argc, char *argv[]) {
        string norm("NONE");
        string unit("BP");
        struct timespec t0,t1,st,en;

        double tol=1.0e-6;
        int maxiter=500;
        int threads = 1;
        int verb = 1;
	string ob("oe");
        int opt;

	clock_gettime(CLOCK_REALTIME,&st);

	while ((opt = getopt(argc, argv, "o:t:I:T:n:v:h")) != -1) {
        	switch (opt) {
                                case 'o':
                                        if (strcmp(optarg,"observed") == 0) ob = "observed";
					else {
						usage(argv[0]);
                                        	exit(EXIT_FAILURE);
					}
                                        break;
                                case 't':
                                        tol = atof(optarg);
                                        break;
                                case 'I':
                                        maxiter=atoi(optarg);
                                        break;
                                case 'n':
                                        norm=optarg;
                                        break;
                                case 'T':
                                        threads=atoi(optarg);
                                        break;
                                case 'v':
                                        verb=atoi(optarg);
                                        break;
                                case 'h':
                                        usage(argv[0]);
                                        exit(EXIT_SUCCESS);
                                default:
                                        usage(argv[0]);
                                        exit(EXIT_FAILURE);
                }
        }

	if (argc - optind < 3) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	string fname = argv[optind++];

	char *out_name = argv[optind++];
	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else if (verb) printf("Writing eigenvectors to %s\n\n", argv[optind-1]);
	int binsize = atoi(argv[optind++]);
	fprintf(fout,"track type=wiggle_0\n");

        hictk::File f(fname, binsize);

        std::vector<int> ii;
	std::vector<int> jj;
	std::vector<double> xx;
	long p,q;
	long m;
	int k;
	int i;

	string hg19("hg19");
	string hg38("hg38");

// chromosome map for finding matrix
    long master = 0L;
    string genomeID;
    int numChromosomes = 0;
    int version = 0;
    long nviPosition = 0;
    long nviLength = 0;
    long totalFileSize;

	if (genomeID != hg19 && genomeID != hg38) if (verb) {
		cout << genomeID;
		cout << " is not currently supported; no sign flip will be attempted!" << endl << endl;
	}

	int iter0 = 0;

//	loop accross all chromosomes
	for (const auto& chrom: f.chromosomes()) {
                const std::string chrom_name{chrom.name()};
		if (chrom_name == "Y" || chrom_name == "M" || chrom_name == "MT") continue;
		if (chrom_name == "chrY" || chrom_name == "chrM" || chrom_name == "chrMT") continue;
		if (chrom.is_all()) continue;
		clock_gettime(CLOCK_REALTIME,&t0);
                ii.clear();
                jj.clear();
                xx.clear();
                const auto sel = f.fetch(chrom.name(), hictk::balancing::Method{norm});
                const auto offset = f.bins().at(chrom, 0).id();
                std::for_each(sel.begin<double>(), sel.end<double>(),
                    [&](const hictk::ThinPixel<double>&p){
                                ii.push_back(p.bin1_id - offset);
                                jj.push_back(p.bin2_id - offset);
                                if (p.bin1_id == p.bin2_id) {
                                  xx.push_back(p.count * 0.5);
                                } else {
                                  xx.push_back(p.count);
                                }

                });

                ii.shrink_to_fit();
                jj.shrink_to_fit();
                xx.shrink_to_fit();

                if (f.is_hic()) {
                  f.get<hictk::hic::File>().clear_cache();
                  f.get<hictk::hic::File>().purge_footer_cache();
                }

                if (ii.empty()) {
                  continue;
                  // exit(EXIT_FAILURE);
                }
		clock_gettime(CLOCK_REALTIME,&t1);
		if (verb > 1) printf("chromosome %s:  took %10.3f seconds for %ld records\n",chrom_name.c_str(),((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9,m);
		m = ii.size();

		int N = (int) ceil(chrom.size()/((double) binsize));

		double l1,l2,er;
		double *a1 = (double *) malloc(N*sizeof(double));
		double *a2 = (double *) malloc(N*sizeof(double));
		double *ev = (double *) malloc(N*sizeof(double));

		long init_seed = 1234563;
		srand48(init_seed);
		for (p=0;p<N;p++) {
		        a1[p] = drand48();
		        a2[p] = drand48();
		}

		clock_gettime(CLOCK_REALTIME,&t0);
		int iter;
		iter = bestEigen(ii.size(),ii.data(),jj.data(),xx.data(),&N,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
		clock_gettime(CLOCK_REALTIME,&t1);
		if (verb > 1) {
			printf("total %d iterations\n",iter);
                	printf("iterations took %10.3f seconds\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9);
			printf("lam1 = %g; lam2 = %g; lam1/lam2 = %g; er = %g; error in EV = %g\n",l1,l2,l1/l2,er,er/(l1-l2));
			printf("                           -------------------- \n");
		}

		fprintf(fout,"fixedStep ");
                char *chr1 = (char *) malloc((strlen(chrom_name.c_str())+4)*sizeof(char));
                if (!strstr(chrom_name.c_str(),"chr")) strcpy(chr1,"chr");
                else strcpy(chr1,"");
                strcat(chr1,chrom_name.c_str());
                if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
                fprintf(fout,"chrom=%s ",chr1);
                fprintf(fout,"start=1 step=%d span=%d\n",binsize,binsize);
		char *genome1 = const_cast<char*> (genomeID.c_str());
		if (100000 % binsize == 0) int junk = flipSign(genome1,ev,N,chr1,binsize);
                for (int j=0;j< N;j++) {
			if (j == N-1) fprintf(fout,"fixedStep chrom=%s start=%d step=%ld span=%ld\n",chr1,j*binsize+1,chrom.size() % binsize, chrom.size() % binsize );
			if (!isnan(ev[j])) fprintf(fout,"%f\n",ev[j]);
			else fprintf(fout,"%s\n","0");
//			else fprintf(fout,"%s\n","NaN");
		}
                fflush(fout);
	}
	fclose(fout);
	clock_gettime(CLOCK_REALTIME,&en);
	if (verb) printf("\n**************    all together took %10.3f seconds\n",((double) (en.tv_sec - st.tv_sec)) + ((double) (en.tv_nsec - st.tv_nsec))/1e9);
	return(0);
}

