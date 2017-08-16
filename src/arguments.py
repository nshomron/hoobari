import argparse
import os

tmp_dir = os.path.join(os.getcwd(), 'tmp_hb')

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = tmp_dir, help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')
parser.add_argument("-r", "--region", default = False, help = "run on a specific region as explained in pyvcf documentation")
parser.add_argument("-v", "--verbosity", action = 'store_true', help = "Prints more detailed debugging information")
parser.add_argument("-l", "--plot_lengths", default = False, action = 'store_true', help = "Creates a plot of the length distributions")
parser.add_argument("-d", "--db", default = os.path.join(tmp_dir, 'hoobari.db'), help = 'db path')
parser.add_argument("-D", "--db_prefix", default = 'False', help = '''If hoobari is run split, all sqlite databases are expected to be in the
									same location as the processed database, sharing some unique prefix''')
parser.add_argument("-model", "--model", default = 'simple', help = '	model for likelihoods calculation. possible values: "simple" \
									(Bayesian model based only on fetal fraction and parental genotypes), \
									"lengths" (use different fetal fraction per fragment length), \
									"origin" (use fragments that are very likely to be fetal, \
									based on other SNPs on these fragments)') # TODO: remove this argument and change the likelihoods function before relsease version

args = parser.parse_args()