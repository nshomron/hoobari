import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--maternal_sample_name", help = 'maternal sample name as appears in parents vcf')
parser.add_argument("-p", "--paternal_sample_name", help = 'paternal sample name as appears in parents vcf')
parser.add_argument("-f", "--fetal_sample_name", help = 'fetal sample name to write in the outputvcf')
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-cfdna_vcf", "--cfdna_vcf", help = 'The maternal plasma cfDNA VCF file')
parser.add_argument("-t", "--tmp_dir", default = os.path.join(os.getcwd(), 'tmp_hb'), help = 'Directory for temporary files')
parser.add_argument("-o", "--vcf_output", default = False, help = 'path for vcf output')
parser.add_argument("-r", "--region", default = False, help = "run on a specific region as explained in pyvcf documentation")
parser.add_argument("-v", "--verbosity", action = 'store_true', help = "Prints more detailed debugging information")
parser.add_argument("-dbt", "--dbtype", default = 'mysql', help = 'mysql or sqlite')
parser.add_argument("-mysql", "--mysql", default = 'tomr@nshomron.tau.ac.il/var/opt/rocks/mysql/mysql.sock:40000', help = 'sqlserver connection information - user@host/path/to/socket.sock:port')
parser.add_argument("-db", "--db", default = 'hoobari', help = 'db name if mySQL or db path if SQLite')
parser.add_argument("-model", "--model", default = 'simple', help = '	model for likelihoods calculation. possible values: "simple" \
									(Bayesian model based only on fetal fraction and parental genotypes), \
									"lengths" (use different fetal fraction per fragment length), \
									"origin" (use fragments that are very likely to be fetal, \
									based on other SNPs on these fragments)') # TODO: remove this argument and change the likelihoods function before relsease version

args = parser.parse_args()