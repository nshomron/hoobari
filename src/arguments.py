import argparse
import os

tmp_dir = os.path.join(os.getcwd(), 'tmp_hb')

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fetal_sample_name", default = 'FETUS', help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = tmp_dir, help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')
parser.add_argument("-r", "--region", default = False, help = "run on a specific region as explained in pyvcf documentation")
parser.add_argument("-v", "--verbosity", action = 'store_true', help = "Prints more detailed debugging information")
parser.add_argument("-w", "--window", default = 3, help = "Window size for lengths")
parser.add_argument("-%", "--fetal_fraction", default = False, help = "specify a known fetal fraction to pass hoobari's calculation of it")
parser.add_argument("-P", "--use_prior_ff_dist", default = False, action = 'store_true', help = "whether to user a prior, known length ratios distribution")
parser.add_argument("-l", "--plot_lengths", default = False, action = 'store_true', help = "Creates a plot of the length distributions")
parser.add_argument("-q", "--qnames", default = False, action = 'store_true', help = "Creates lists of fetal and shared qnames")
parser.add_argument("-d", "--db", default = os.path.join(tmp_dir, 'hoobari.db'), help = 'db path')
parser.add_argument("-D", "--db_prefix", default = False, help = '''If hoobari is run split, all sqlite databases are expected to be in the
									same location as the processed database, sharing some unique prefix''')
parser.add_argument("-@", "--cores", default = 1, help = 'number of cores to run pre-processing when run split')
parser.add_argument("-model", "--model", default = 'lengths', help = '	model for likelihoods calculation. possible values: "simple" \
									(Bayesian model based only on fetal fraction and parental genotypes), \
									"lengths" (use different fetal fraction per fragment length), \
									"origin" (use fragments that are very likely to be fetal, \
									based on other SNPs on these fragments)') # TODO: remove this argument and change the likelihoods function before relsease version

args = parser.parse_args()
