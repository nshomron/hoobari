'''
receives err output pf freebayes -dd as stdin, and for every assessed position, prints out the reads supporting each genotype
'''

import sys
import os
import re
import math
import pysam
import vcf
import db
import argparse
import pandas as pd
import numpy as np
import random
import string

# --------- parse args ---------
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam_file", help = 'The maternal cfDNA bam file\'s path')
parser.add_argument("-t", "--tmp_dir", default = 'tmp_hb')
parser.add_argument("-r", "--region", default = False)
parser.add_argument("-s", "--downsample", default = False, help = 'in the format of current_fraction,new_fraction')
parser.add_argument("-d", "--debug", action = 'store_true', default = False, help = """By default, Freebayes' stderr is used by this patch,
						which can cause a problem when actually trying to debug.
						This flag causes this information to be printed. Please
						note that the data, if saved, ends up in a very large file.""")
parser.add_argument("-parents_vcf", "--parents_vcf", help = 'bgzipped vcf of parents, indexed by tabix')
parser.add_argument("-m", "--m_bam", help = 'maternal bam file')
parser.add_argument("-p", "--p_bam", help = 'paternal bam file')
parser.add_argument("-db", "--db", default = 'hoobari', help = 'db name, or db prefix if hoobari is run per region')
parser.add_argument("-origin", "--origin", action = 'store_true', help = 'use this if hoobari will be run with -model origin')
args = parser.parse_args()
# ------------------------------

var_type_dic = {'.': 0, 'snp': 1, 'mnp': 2, 'ins': 3, 'del': 4, 'complex': 5}

if args.downsample:
	# Family G1 length distributions, used for downsampling experiments.
	# If testing on another family, these should be changed accordingly
	maternal_length_distribution = [932420.0, 0.0, 42.0, 30.0, 33.0, 13.0, 5.0, 12.0, 1.0, 10.0, 12.0, 24.0, 5.0, 19.0, 19.0, 11.0, 4.0, 12.0, 10.0, 33.0, 20.0, 24.0, 20.0, 17.0, 20.0, 27.0, 23.0, 0.0, 18.0, 33.0, 202.0, 199.0, 235.0, 289.0, 344.0, 579.0, 595.0, 769.0, 861.0, 999.0, 1175.0, 1259.0, 1320.0, 1494.0, 1741.0, 1844.0, 2210.0, 2560.0, 3007.0, 3463.0, 3944.0, 4086.0, 4192.0, 4252.0, 4690.0, 5077.0, 6127.0, 6605.0, 7886.0, 9226.0, 9692.0, 10163.0, 10194.0, 10131.0, 10595.0, 11040.0, 12257.0, 13148.0, 14418.0, 16167.0, 17648.0, 15128.0, 15007.0, 16971.0, 17467.0, 18373.0, 19151.0, 20139.0, 21177.0, 26164.0, 29198.0, 28791.0, 26212.0, 23261.0, 22058.0, 21641.0, 22323.0, 21955.0, 23707.0, 27277.0, 32932.0, 35349.0, 34879.0, 29850.0, 25276.0, 23343.0, 23827.0, 23700.0, 26830.0, 28956.0, 34220.0, 20898.0, 36275.0, 35450.0, 33648.0, 31283.0, 30862.0, 27866.0, 35024.0, 40057.0, 47836.0, 56644.0, 54229.0, 48657.0, 44156.0, 44137.0, 45161.0, 47678.0, 52728.0, 61226.0, 76802.0, 98109.0, 116008.0, 111493.0, 89837.0, 73537.0, 69108.0, 71856.0, 78988.0, 88842.0, 106582.0, 135258.0, 177807.0, 219904.0, 227725.0, 203615.0, 181119.0, 185468.0, 215160.0, 265173.0, 330380.0, 389963.0, 420044.0, 427499.0, 443835.0, 461171.0, 465263.0, 475301.0, 502465.0, 544459.0, 590012.0, 655360.0, 703573.0, 711397.0, 698999.0, 702389.0, 744736.0, 816905.0, 908628.0, 1006799.0, 1088799.0, 1149919.0, 1211761.0, 1318346.0, 1477111.0, 1638354.0, 1725063.0, 1762417.0, 1758499.0, 1728113.0, 1673014.0, 1609882.0, 1565129.0, 1559369.0, 1578846.0, 1587192.0, 1576801.0, 1517744.0, 1425050.0, 1320479.0, 1212346.0, 1117479.0, 1037282.0, 976727.0, 932990.0, 903709.0, 879053.0, 851054.0, 810485.0, 759722.0, 699627.0, 637289.0, 583142.0, 534483.0, 497141.0, 468319.0, 444912.0, 424123.0, 408201.0, 384643.0, 360415.0, 332749.0, 305352.0, 278950.0, 256441.0, 236090.0, 221691.0, 208086.0, 196031.0, 184458.0, 172171.0, 159028.0, 146250.0, 133753.0, 121426.0, 111925.0, 103122.0, 95377.0, 89086.0, 83903.0, 78387.0, 72458.0, 66258.0, 60368.0, 55747.0, 50816.0, 45999.0, 43160.0, 39724.0, 36467.0, 34658.0, 31532.0, 29004.0, 26798.0, 24917.0, 22413.0, 20426.0, 19424.0, 17991.0, 16766.0, 16143.0, 14684.0, 13925.0, 13013.0, 12074.0, 11126.0, 10531.0, 9917.0, 9386.0, 8925.0, 8612.0, 8333.0, 7892.0, 7671.0, 7280.0, 7201.0, 6755.0, 6544.0, 6392.0, 6333.0, 6076.0, 6106.0, 6056.0, 6005.0, 5730.0, 5890.0, 5453.0, 5422.0, 5419.0, 5801.0, 5903.0, 5975.0, 6048.0, 6592.0, 6580.0, 6131.0, 6239.0, 6107.0, 6676.0, 6952.0, 6965.0, 7299.0, 7717.0, 8086.0, 8340.0, 8681.0, 8997.0, 9296.0, 10004.0, 11037.0, 11226.0, 11931.0, 13006.0, 13681.0, 14475.0, 15474.0, 16523.0, 17781.0, 18729.0, 19951.0, 20956.0, 22180.0, 23129.0, 24703.0, 26540.0, 27984.0, 30279.0, 32072.0, 34391.0, 35566.0, 37144.0, 38989.0, 40149.0, 42065.0, 43974.0, 45789.0, 47624.0, 49274.0, 50345.0, 51995.0, 53564.0, 54338.0, 55368.0, 57050.0, 57754.0, 58515.0, 58902.0, 59959.0, 60104.0, 61145.0, 61850.0, 63624.0, 64761.0, 64773.0, 64295.0, 64375.0, 64340.0, 63349.0, 63773.0, 63602.0, 63761.0, 64294.0, 64407.0, 65073.0, 65229.0, 64353.0, 64427.0, 64138.0, 63769.0, 63079.0, 62738.0, 62137.0, 61556.0, 60821.0, 60502.0, 59970.0, 59543.0, 58888.0, 58129.0, 57436.0, 56628.0, 55043.0, 55213.0, 53862.0, 52135.0, 51946.0, 50601.0, 49683.0, 48687.0, 47959.0, 47464.0, 45779.0, 44303.0, 43241.0, 41642.0, 40484.0, 38979.0, 38464.0, 36946.0, 36026.0, 34876.0, 33621.0, 32167.0, 31301.0, 29716.0, 28478.0, 27282.0, 26257.0, 25589.0, 24160.0, 23478.0, 22210.0, 21250.0, 20684.0, 19087.0, 18461.0, 16997.0, 16406.0, 15630.0, 14769.0, 14360.0, 13691.0, 13040.0, 12343.0, 11695.0, 11043.0, 10341.0, 9642.0, 9352.0, 8655.0, 8304.0, 7963.0, 7540.0, 7163.0, 6550.0, 6407.0, 6076.0, 5698.0, 5391.0, 5028.0, 4728.0, 4595.0, 4439.0, 4165.0, 3886.0, 3750.0, 3556.0, 3392.0, 3216.0, 3108.0, 3002.0, 2847.0, 2845.0, 2847.0, 2653.0, 2619.0, 2531.0, 2429.0, 2336.0, 2374.0, 2278.0, 2269.0, 2303.0, 2265.0, 2275.0, 2209.0, 2335.0, 2265.0, 2304.0, 2371.0, 2371.0, 2343.0, 2481.0, 2450.0, 2569.0, 2601.0, 2697.0, 2583.0, 2823.0, 2752.0, 2789.0, 2947.0, 2921.0, 3170.0, 3190.0, 3233.0, 3322.0, 3238.0, 3284.0, 3349.0, 3485.0, 3479.0, 3738.0, 3733.0, 3848.0, 3829.0, 3922.0, 3920.0, 3948.0, 4077.0, 4065.0, 4140.0, 4074.0, 4103.0, 4417.0, 4262.0, 4283.0, 4283.0, 4294.0, 4366.0, 4389.0, 4419.0, 4486.0, 4509.0, 4405.0, 4355.0, 4510.0, 4610.0, 4545.0, 4626.0, 4561.0, 4747.0, 4443.0, 4691.0, 4588.0, 4510.0, 4715.0, 4745.0, 4758.0, 4822.0, 5005.0, 4830.0, 4755.0, 4860.0, 4979.0, 4885.0, 5033.0, 5021.0, 5088.0, 5021.0, 5197.0, 5228.0, 5088.0, 5139.0, 5156.0, 5144.0, 5176.0, 5218.0, 5224.0, 5187.0, 5378.0, 5372.0, 5340.0, 5174.0, 5132.0, 5177.0, 5134.0, 5151.0, 5225.0, 5152.0, 5070.0, 5050.0, 4999.0, 4980.0, 5014.0, 4843.0, 4918.0, 4999.0, 4801.0, 4861.0, 4885.0, 4765.0, 4679.0, 4664.0, 4529.0, 4486.0, 4525.0, 4392.0, 4214.0, 4245.0, 4025.0, 4089.0, 3980.0, 3946.0, 3931.0, 3824.0, 3806.0, 3603.0, 3539.0, 3392.0, 3306.0, 3275.0, 3239.0, 3185.0, 3182.0, 3034.0, 2924.0, 2888.0, 2845.0, 2698.0, 2627.0, 2562.0, 2378.0, 2400.0, 2350.0, 2256.0, 2128.0, 2134.0, 2058.0, 2049.0, 1907.0, 1927.0, 1885.0, 1818.0, 1732.0, 1639.0, 1597.0, 1515.0, 1477.0, 1494.0, 1322.0, 1288.0, 1304.0, 1274.0, 1210.0, 1171.0, 1077.0, 1148.0, 1051.0, 992.0, 1001.0, 968.0, 1033.0, 928.0, 970.0, 871.0, 853.0, 836.0, 788.0, 787.0, 836.0, 770.0, 749.0, 760.0, 760.0, 792.0, 734.0, 694.0, 751.0, 714.0, 669.0, 681.0, 714.0, 724.0, 679.0, 656.0, 698.0, 712.0, 687.0, 674.0, 651.0, 645.0, 672.0, 617.0, 645.0, 663.0, 630.0, 624.0, 685.0, 689.0, 643.0, 632.0, 661.0, 636.0, 654.0, 604.0, 670.0, 627.0, 678.0, 640.0, 624.0, 640.0, 675.0, 632.0, 598.0, 611.0, 598.0, 632.0, 630.0, 605.0, 677.0, 681.0, 617.0, 638.0, 638.0, 625.0, 657.0, 614.0, 682.0, 638.0, 625.0, 639.0, 634.0, 706.0, 654.0, 676.0, 709.0, 662.0, 692.0, 715.0, 656.0, 682.0, 705.0, 654.0, 748.0, 679.0, 739.0, 714.0, 771.0, 684.0, 697.0, 703.0, 741.0, 788.0, 691.0, 705.0, 746.0, 743.0, 768.0, 743.0, 718.0, 734.0, 771.0, 718.0, 736.0, 719.0, 737.0, 721.0, 744.0, 785.0, 789.0, 725.0, 739.0, 723.0, 777.0, 794.0, 733.0, 778.0, 736.0, 686.0, 754.0, 743.0, 726.0, 709.0, 699.0, 747.0, 746.0, 759.0, 696.0, 672.0, 707.0, 680.0, 663.0, 651.0, 729.0, 671.0, 685.0, 643.0, 599.0, 629.0, 542.0, 624.0, 629.0, 623.0, 634.0, 568.0, 605.0, 550.0, 554.0, 559.0, 538.0, 533.0, 534.0, 528.0, 487.0, 488.0, 485.0, 497.0, 464.0, 482.0, 465.0, 421.0, 460.0, 397.0, 391.0, 422.0, 425.0, 438.0, 429.0, 335.0, 352.0, 378.0, 346.0, 322.0, 343.0, 332.0, 365.0, 330.0, 320.0, 300.0, 320.0, 303.0, 292.0, 268.0, 256.0, 262.0, 250.0, 262.0, 243.0, 251.0, 208.0, 233.0, 205.0, 232.0, 195.0, 193.0, 187.0, 206.0, 186.0, 186.0, 195.0, 184.0, 190.0, 191.0, 201.0, 190.0, 0.0, 156.0, 164.0, 138.0, 173.0, 177.0, 152.0, 129.0, 152.0, 163.0, 145.0, 144.0, 133.0, 149.0, 137.0, 145.0, 129.0, 133.0, 135.0, 137.0, 145.0, 147.0, 115.0, 116.0, 147.0, 153.0, 145.0, 146.0, 131.0, 138.0, 123.0, 133.0, 129.0, 138.0, 126.0, 108.0, 129.0, 125.0, 118.0, 143.0, 126.0, 121.0, 127.0, 118.0, 122.0, 139.0, 124.0, 109.0, 117.0, 136.0, 119.0, 131.0, 122.0, 115.0, 122.0, 129.0, 120.0, 103.0, 133.0, 114.0, 143.0, 133.0, 153.0, 128.0, 126.0, 129.0, 103.0, 156.0, 123.0, 152.0, 118.0, 128.0, 123.0, 95.0, 125.0, 130.0, 124.0, 149.0, 144.0, 135.0, 143.0, 136.0, 155.0, 113.0, 132.0, 120.0, 138.0, 122.0, 131.0, 0.0, 144.0, 127.0, 124.0, 144.0, 129.0, 139.0, 130.0, 153.0, 130.0, 156.0, 126.0, 145.0, 120.0, 144.0, 123.0, 144.0, 126.0, 132.0, 120.0, 139.0, 130.0, 143.0, 138.0, 109.0, 143.0, 120.0, 117.0, 145.0, 110.0, 115.0, 116.0, 113.0, 124.0, 129.0, 112.0, 115.0, 95.0, 131.0, 97.0, 98.0, 0.0, 112.0, 104.0, 104.0, 111.0, 101.0, 91.0, 92.0, 100.0, 104.0, 103.0, 92.0, 95.0, 103.0, 87.0, 84.0, 85.0, 81.0, 79.0, 84.0, 0.0, 0.0, 74.0, 85.0, 64.0, 74.0, 78.0, 0.0, 73.0, 53.0, 0.0, 57.0, 60.0, 70.0, 41.0, 64.0, 59.0, 55.0, 62.0, 56.0, 50.0, 46.0, 61.0, 0.0, 56.0, 53.0, 48.0, 54.0, 51.0, 44.0, 40.0, 44.0, 42.0, 39.0, 53.0]
	fetal_length_distribution = [263692.0, 0.0, 21.0, 12.0, 18.0, 16.0, 12.0, 7.0, 10.0, 8.0, 5.0, 5.0, 7.0, 4.0, 9.0, 5.0, 5.0, 8.0, 12.0, 6.0, 11.0, 6.0, 7.0, 7.0, 11.0, 13.0, 10.0, 37.0, 19.0, 17.0, 67.0, 84.0, 89.0, 87.0, 121.0, 189.0, 321.0, 291.0, 339.0, 421.0, 433.0, 411.0, 475.0, 488.0, 553.0, 767.0, 770.0, 957.0, 1142.0, 1392.0, 1466.0, 1451.0, 1373.0, 1420.0, 1504.0, 1818.0, 2119.0, 2854.0, 3477.0, 4010.0, 4353.0, 4333.0, 4191.0, 4198.0, 4601.0, 5096.0, 5055.0, 5604.0, 6726.0, 8177.0, 9290.0, 10699.0, 10012.0, 8044.0, 8021.0, 8526.0, 8862.0, 9821.0, 12559.0, 15108.0, 18505.0, 18953.0, 15393.0, 12300.0, 11121.0, 11231.0, 11771.0, 13051.0, 14413.0, 18289.0, 24523.0, 29741.0, 29345.0, 24697.0, 20937.0, 18995.0, 17655.0, 18529.0, 20290.0, 26573.0, 31318.0, 50749.0, 35090.0, 32228.0, 29037.0, 26758.0, 26934.0, 33307.0, 31643.0, 38481.0, 48613.0, 58844.0, 60077.0, 53335.0, 46140.0, 44425.0, 46107.0, 48404.0, 52334.0, 61747.0, 79937.0, 104344.0, 119780.0, 115354.0, 93484.0, 78834.0, 76633.0, 79983.0, 84404.0, 91789.0, 107276.0, 131108.0, 163364.0, 192969.0, 195081.0, 171526.0, 148700.0, 142080.0, 148397.0, 168324.0, 200130.0, 223601.0, 236131.0, 237978.0, 238973.0, 236299.0, 225479.0, 211707.0, 205471.0, 204055.0, 214931.0, 225050.0, 237476.0, 237569.0, 227774.0, 220481.0, 218767.0, 218950.0, 224067.0, 230383.0, 238416.0, 243665.0, 243998.0, 247023.0, 252927.0, 257365.0, 257149.0, 250101.0, 245670.0, 236804.0, 227532.0, 216921.0, 208696.0, 200484.0, 194565.0, 190116.0, 183568.0, 174387.0, 164173.0, 153013.0, 141889.0, 131867.0, 122884.0, 115228.0, 108483.0, 103148.0, 97632.0, 92300.0, 86967.0, 81052.0, 74917.0, 68876.0, 63679.0, 58514.0, 54762.0, 50798.0, 48619.0, 46111.0, 42906.0, 40254.0, 37197.0, 34628.0, 31957.0, 29794.0, 27690.0, 25698.0, 24132.0, 22467.0, 21362.0, 19701.0, 18154.0, 16863.0, 15878.0, 14296.0, 13094.0, 12091.0, 11481.0, 10919.0, 9885.0, 9301.0, 8818.0, 8171.0, 7717.0, 7138.0, 6562.0, 5711.0, 5577.0, 4992.0, 4809.0, 4361.0, 4257.0, 4030.0, 3972.0, 3727.0, 3193.0, 2969.0, 2907.0, 2534.0, 2516.0, 2501.0, 2553.0, 2605.0, 2496.0, 2367.0, 2160.0, 1978.0, 2010.0, 1875.0, 1832.0, 1914.0, 2001.0, 2148.0, 2274.0, 2281.0, 2214.0, 1933.0, 1949.0, 1815.0, 1895.0, 2023.0, 2333.0, 2543.0, 2685.0, 2920.0, 2956.0, 2883.0, 2762.0, 2599.0, 2535.0, 2605.0, 2759.0, 3185.0, 3603.0, 4373.0, 4382.0, 4363.0, 3507.0, 3407.0, 3380.0, 3918.0, 3514.0, 3476.0, 3906.0, 4132.0, 4118.0, 4185.0, 4197.0, 4174.0, 4231.0, 4297.0, 4572.0, 4825.0, 5102.0, 5158.0, 5394.0, 5281.0, 5242.0, 5417.0, 5442.0, 5825.0, 6155.0, 6268.0, 6437.0, 6400.0, 6419.0, 6478.0, 6794.0, 6774.0, 6908.0, 7193.0, 7185.0, 7231.0, 7225.0, 7303.0, 7396.0, 7486.0, 7297.0, 7436.0, 7393.0, 7240.0, 7141.0, 7019.0, 7226.0, 6985.0, 6926.0, 6804.0, 6724.0, 6559.0, 6431.0, 6287.0, 6191.0, 6136.0, 5878.0, 5779.0, 5698.0, 5611.0, 5309.0, 5249.0, 5163.0, 4995.0, 4845.0, 4660.0, 4657.0, 4478.0, 4377.0, 4185.0, 3972.0, 4028.0, 3893.0, 3711.0, 3678.0, 3621.0, 3500.0, 3336.0, 3284.0, 3181.0, 3128.0, 2933.0, 2895.0, 2883.0, 2695.0, 2574.0, 2586.0, 2457.0, 2330.0, 2292.0, 2158.0, 2129.0, 2126.0, 2028.0, 1894.0, 1876.0, 1714.0, 1679.0, 1577.0, 1583.0, 1558.0, 1425.0, 1410.0, 1388.0, 1338.0, 1341.0, 1249.0, 1161.0, 1081.0, 1033.0, 1023.0, 975.0, 914.0, 911.0, 811.0, 828.0, 772.0, 750.0, 673.0, 698.0, 653.0, 632.0, 591.0, 566.0, 537.0, 530.0, 473.0, 517.0, 460.0, 408.0, 381.0, 376.0, 393.0, 382.0, 365.0, 359.0, 318.0, 334.0, 301.0, 311.0, 282.0, 291.0, 277.0, 306.0, 262.0, 251.0, 248.0, 230.0, 243.0, 256.0, 250.0, 283.0, 237.0, 255.0, 247.0, 258.0, 260.0, 236.0, 278.0, 259.0, 285.0, 258.0, 233.0, 243.0, 260.0, 259.0, 251.0, 283.0, 270.0, 298.0, 279.0, 277.0, 326.0, 317.0, 305.0, 331.0, 316.0, 324.0, 312.0, 313.0, 329.0, 340.0, 317.0, 335.0, 342.0, 343.0, 344.0, 331.0, 385.0, 337.0, 321.0, 335.0, 352.0, 345.0, 340.0, 337.0, 354.0, 318.0, 340.0, 327.0, 331.0, 335.0, 338.0, 328.0, 316.0, 316.0, 323.0, 343.0, 322.0, 318.0, 299.0, 310.0, 292.0, 289.0, 291.0, 326.0, 294.0, 276.0, 298.0, 298.0, 278.0, 249.0, 250.0, 305.0, 271.0, 263.0, 276.0, 263.0, 259.0, 231.0, 225.0, 247.0, 243.0, 224.0, 249.0, 229.0, 222.0, 220.0, 258.0, 214.0, 212.0, 205.0, 205.0, 211.0, 208.0, 191.0, 186.0, 214.0, 190.0, 190.0, 172.0, 177.0, 194.0, 193.0, 205.0, 188.0, 183.0, 175.0, 194.0, 186.0, 174.0, 178.0, 173.0, 155.0, 170.0, 137.0, 146.0, 177.0, 176.0, 158.0, 145.0, 145.0, 155.0, 125.0, 130.0, 132.0, 138.0, 112.0, 120.0, 123.0, 102.0, 137.0, 122.0, 137.0, 107.0, 104.0, 104.0, 99.0, 93.0, 93.0, 98.0, 91.0, 100.0, 78.0, 85.0, 86.0, 84.0, 91.0, 89.0, 60.0, 72.0, 72.0, 69.0, 81.0, 77.0, 62.0, 79.0, 60.0, 55.0, 57.0, 51.0, 54.0, 65.0, 59.0, 68.0, 61.0, 38.0, 56.0, 41.0, 60.0, 46.0, 51.0, 59.0, 35.0, 52.0, 51.0, 33.0, 47.0, 52.0, 49.0, 45.0, 28.0, 43.0, 44.0, 46.0, 37.0, 35.0, 39.0, 40.0, 42.0, 33.0, 33.0, 37.0, 43.0, 32.0, 35.0, 45.0, 33.0, 31.0, 32.0, 24.0, 35.0, 38.0, 28.0, 37.0, 39.0, 37.0, 25.0, 45.0, 40.0, 33.0, 33.0, 35.0, 36.0, 35.0, 25.0, 28.0, 36.0, 36.0, 21.0, 40.0, 28.0, 28.0, 30.0, 26.0, 29.0, 32.0, 25.0, 34.0, 34.0, 28.0, 33.0, 15.0, 33.0, 26.0, 31.0, 31.0, 35.0, 34.0, 26.0, 30.0, 19.0, 23.0, 22.0, 17.0, 18.0, 19.0, 21.0, 23.0, 28.0, 31.0, 26.0, 25.0, 21.0, 39.0, 30.0, 21.0, 23.0, 26.0, 18.0, 27.0, 25.0, 22.0, 18.0, 27.0, 19.0, 20.0, 23.0, 18.0, 23.0, 21.0, 25.0, 20.0, 27.0, 24.0, 19.0, 17.0, 15.0, 21.0, 25.0, 23.0, 22.0, 19.0, 22.0, 20.0, 15.0, 28.0, 27.0, 9.0, 20.0, 30.0, 26.0, 15.0, 17.0, 16.0, 23.0, 14.0, 16.0, 18.0, 17.0, 14.0, 21.0, 18.0, 16.0, 11.0, 19.0, 12.0, 11.0, 18.0, 8.0, 17.0, 16.0, 17.0, 21.0, 16.0, 17.0, 12.0, 15.0, 11.0, 9.0, 16.0, 12.0, 13.0, 16.0, 15.0, 12.0, 16.0, 19.0, 10.0, 14.0, 20.0, 9.0, 12.0, 11.0, 9.0, 12.0, 12.0, 14.0, 17.0, 14.0, 15.0, 7.0, 11.0, 14.0, 7.0, 10.0, 8.0, 10.0, 16.0, 10.0, 9.0, 11.0, 8.0, 5.0, 9.0, 10.0, 10.0, 5.0, 14.0, 9.0, 7.0, 9.0, 7.0, 6.0, 8.0, 10.0, 4.0, 5.0, 5.0, 10.0, 6.0, 8.0, 9.0, 5.0, 9.0, 3.0, 8.0, 5.0, 9.0, 7.0, 4.0, 7.0, 2.0, 4.0, 4.0, 4.0, 4.0, 3.0, 5.0, 6.0, 0.0, 5.0, 5.0, 9.0, 7.0, 4.0, 5.0, 6.0, 5.0, 4.0, 4.0, 2.0, 3.0, 2.0, 6.0, 3.0, 6.0, 6.0, 6.0, 6.0, 3.0, 3.0, 6.0, 5.0, 1.0, 6.0, 9.0, 7.0, 3.0, 5.0, 3.0, 7.0, 3.0, 2.0, 10.0, 10.0, 2.0, 3.0, 3.0, 3.0, 3.0, 7.0, 4.0, 7.0, 5.0, 1.0, 1.0, 3.0, 4.0, 1.0, 7.0, 1.0, 5.0, 2.0, 6.0, 6.0, 5.0, 4.0, 2.0, 3.0, 6.0, 6.0, 2.0, 2.0, 4.0, 1.0, 6.0, 5.0, 3.0, 1.0, 3.0, 4.0, 4.0, 7.0, 7.0, 8.0, 2.0, 4.0, 6.0, 5.0, 6.0, 1.0, 3.0, 6.0, 5.0, 4.0, 5.0, 4.0, 5.0, 0.0, 3.0, 6.0, 3.0, 1.0, 3.0, 3.0, 2.0, 4.0, 4.0, 3.0, 2.0, 1.0, 3.0, 3.0, 3.0, 2.0, 4.0, 2.0, 3.0, 3.0, 4.0, 1.0, 2.0, 2.0, 1.0, 2.0, 5.0, 3.0, 2.0, 4.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0, 9.0, 6.0, 4.0, 0.0, 2.0, 4.0, 2.0, 5.0, 2.0, 5.0, 3.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 4.0, 3.0, 2.0, 4.0, 4.0, 0.0, 0.0, 2.0, 4.0, 1.0, 4.0, 2.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 3.0, 1.0, 1.0, 2.0, 1.0, 4.0, 2.0, 1.0, 2.0, 3.0, 0.0, 2.0, 3.0, 3.0, 4.0, 2.0, 1.0, 3.0, 1.0, 1.0, 3.0, 5.0]
	fetal_frequencies = list(fetal_length_distribution / np.sum(fetal_length_distribution))
	maternal_frequencies = list(maternal_length_distribution / np.sum(maternal_length_distribution))

# --------- functions ---------
def get_parental_genotypes(parents_reader, parental_samples, chrom, position):
	maternal_sample_name, paternal_sample_name = parental_samples
	n_rec = 0
	for rec in parents_reader.fetch(chrom, int(position) - 1, int(position)):
		maternal_gt = rec.genotype(maternal_sample_name).data.GT
		paternal_gt = rec.genotype(paternal_sample_name).data.GT
		n_rec += 1
		if n_rec > 1:
			sys.exit('more than one parental variant in the same position')
	if n_rec == 0:
		maternal_gt, paternal_gt = None, None
	return maternal_gt, paternal_gt

def get_reads_tlen(bam_reader, chrom, position):
	start_pos = max(0, int(position) - 1000)
	end_pos = min(bam_reader.lengths[bam_reader.references.index(chrom)], int(position) + 1000)
	bam_records_at_position = bam_reader.fetch(chrom, start_pos, end_pos) # include a flanking region, since there's local realignment
	tlen_at_position_dic = {}
	for rec in bam_records_at_position:
		tlen_at_position_dic[rec.query_name] = math.fabs(int(rec.template_length))
	return tlen_at_position_dic

def get_fetal_allele_type(maternal_gt, paternal_gt):
	if maternal_gt == '0/0' and paternal_gt in ('0/1','1/1'):
		return 'alt'
	elif maternal_gt == '1/1' and paternal_gt in ('0/0','0/1'):
		return 'ref'
	else:
		return False

def is_fetal_fragment(genotype, ref, alt, fetal_allele = False):

	if ((genotype == ref) and fetal_allele == 'ref') or ((genotype == alt) and fetal_allele == 'alt'):
		return 1
	elif ((genotype == alt) and fetal_allele == 'ref') or ((genotype == ref) and fetal_allele == 'alt'):
		return 0
	else:
		return None

def get_parental_samples_names(m_bam, p_bam):
	parents_sample_names = []
	for bam_file in (m_bam, p_bam):
		bam_file_reader = pysam.AlignmentFile(os.path.join(bam_file), 'rb')
		sample_name = bam_file_reader.header['RG'][0]['SM']
		parents_sample_names.append(sample_name)
		bam_file_reader.close()
	return parents_sample_names

def use_for_fetal_fraction_calculation(maternal_gt, paternal_gt, var_type, is_fetal):
	var_in_ff_positions = (maternal_gt == '0/0' and paternal_gt == '1/1') or (maternal_gt == '1/1' and paternal_gt == '0/0')
	var_is_snp = var_type == 1

	if var_in_ff_positions and var_is_snp:
		if is_fetal == 1:
			return 1
		elif is_fetal == 0:
			return 2
		else:
			return 0
	else:
		return 0 # not for ff

def which_allele_is_fetal(maternal_gt, paternal_gt, ref, alt):
	if maternal_gt == '0/0' and paternal_gt in ('0/1','1/1'):
		return alt
	elif maternal_gt == '1/1' and paternal_gt in ('0/0','0/1'):
		return ref
	else:
		return 'ambiguous'

def qname_generator(size=6, chars=string.ascii_uppercase + string.digits):
	return(''.join(random.choice(chars) for _ in range(size)))

def random_round(x):
	floor_or_ceiling = random.randint(0,1)
	if floor_or_ceiling == 0:
		return(int(np.floor(x)))
	else:
		return(int(np.ceil(x)))

def downsample_position(position_list,
			F,
			D,
			fetal_frequencies,
			maternal_frequencies,
			ref,
			alt,
			fetal_allele):
	'''
	N = original number of fetal reads 
	n = (D/F) * reads with fetal allele
	sample n reads from reads with fetal allele
	for (N - n) reads with shortest isize from the shared list, remove isize and give it one from the maternal size distribution
	add (N - n) reads to the shared allele:
	 - qnames are generated with function
	 - isizes from the distribution of (shared lengths - fetal lengths)
	'''

	if fetal_allele == 'ambiguous':
		# drow n_new_isizes reads from the total allele list using the fetal length
		# distribution, and replace their isizes with isizes from the maternal length distribution
		
		n_new_isizes = random_round((D/F) * len(position_list))
		
		if n_new_isizes > 0:

			position_list_fetal_frequencies = []
			for l in position_list:
				isize = l[1]
				# print(isize, file=sys.stderr)
				if isize >= len(fetal_frequencies):
					isize = np.random.choice(1001, 1, p = fetal_frequencies)[0]
				# print(isize, file=sys.stderr)
				# print(len(fetal_frequencies), file=sys.stderr)
				position_list_fetal_frequencies.append(fetal_frequencies[int(isize)])
			position_list_fetal_frequencies /= np.sum(position_list_fetal_frequencies) # normalize

			new_isizes_idxs = np.random.choice(	len(position_list_fetal_frequencies),
								size = n_new_isizes,
								replace = False,
								p = position_list_fetal_frequencies)
			# print(new_isizes_idxs, file=sys.stderr)
			new_isizes = np.random.choice(1001, size = n_new_isizes, p = maternal_frequencies)
			
			j = 0
			for i in range(len(position_list)):
				if i in new_isizes_idxs:
					position_list[i][1] = new_isizes[j]
					j += 1

		return(position_list)

	else:
		
		# split to two lists
		shared_position_list = []
		fetal_position_list = []
		for line in position_list:
			if line[0] == fetal_allele:
				fetal_position_list.append(line)
			else:
				shared_position_list.append(line)
		
		# assemble fetal list
		n_fetal_original = len(fetal_position_list)
		F_at_site = (2 * n_fetal_original) / len(position_list)
		if F_at_site < 1 and F_at_site > D:
			F = F_at_site
		n_fetal_new = random_round((D/F) * n_fetal_original)
		n_swap = n_fetal_original - n_fetal_new

		fetal_final_list = random.sample(fetal_position_list, n_fetal_new)

		# drow n_swap reads from the shared allele list using the fetal length
		# distribution, and replace their isizes with isizes from the maternal length distribution
		n_new_isizes = min(len(shared_position_list), n_swap)

		if n_new_isizes > 0:

			shared_list_fetal_frequencies = []
			for l in shared_position_list:
				isize = l[1]
				if isize >= len(fetal_frequencies):
					isize = np.random.choice(1001, 1, p = fetal_frequencies)[0]
				shared_list_fetal_frequencies.append(fetal_frequencies[int(isize)])
			shared_list_fetal_frequencies /= np.sum(shared_list_fetal_frequencies) # normalize

			new_isizes_idxs = np.random.choice(	len(shared_list_fetal_frequencies),
								size = n_new_isizes,
								replace = False,
								p = shared_list_fetal_frequencies)

			new_isizes = np.random.choice(1001, size = n_new_isizes, p = maternal_frequencies)
			
			j = 0
			for i in range(len(shared_position_list)):
				if i in new_isizes_idxs:
					shared_position_list[i][1] = new_isizes[j]
					j += 1
		
		# create new records for the maternal list
		# new alleles
		maternal_allele = [i for i in (ref, alt) if i != fetal_allele]
		genos = maternal_allele * n_swap
		# new isizes
		isizes = np.random.choice(1001, size = n_swap, p = maternal_frequencies)
		# new qnames
		qnames = [qname_generator(37) for i in range(n_swap)]

		for geno, isize, qname in zip(genos, isizes, qnames):
			shared_position_list.append([geno, float(isize), qname])
	
		return(shared_position_list + fetal_final_list)		


'''
This patch uses freebayes' algorithm to create a folders' tree: tmp_folder/jsons/chr[1-22,X,Y].
In each chromosome's folder it creates json files named after the position of the variant they represent.

Usage (in bash, using Python 3.5):
freebayes -dd -f [REFERENCE] [OPTIONS] [cfDNA_BAM_FILE] 2>&1 >[OUTPUT] | python freebayes_dd_hoobari_patch.py -b cfDNA_BAM_FILE -t tmp_folder

Explanation:
1) freebayes has to be run with -dd flag, which print more verbose debugging output (and requires "make DEBUG" for installation)
2) Since freebayes' debug information is written to stderr, 2>&1 redirects it to stdout in order to pipe it
3) the tmp_folder created here is the same one you should later use when you run hoobari

1 - var_type - 5 + null
2 - is_fetal - 0/1
3 - ff - 0, 1, 2 (2 + null)

'''

# Initiate variants database
bam_reader = pysam.AlignmentFile(os.path.join(args.bam_file), 'rb')

parents_reader = vcf.Reader(filename = args.parents_vcf)
if args.region:
	dbpath = os.path.join(args.tmp_dir, str(args.region) + '.db')
else:
	dbpath = os.path.join(args.tmp_dir, args.db + '.db')
vardb = db.Variants(dbpath = dbpath)

# get parental sample names
parents_sample_names = get_parental_samples_names(args.m_bam, args.p_bam)

# create sample table in the db
vardb.create_samples_table(parents_sample_names)

# write_variant = False
for line in sys.stdin:
	
	if args.debug:
		print(line, file = sys.stderr, end = '')

	
	if line.startswith('#'):
		print(line, end = '')

	elif line.startswith('position: '):
		initiate_var = True
		line = line.split()
		var = line[1]

	elif line.startswith('haplo_obs'):
		if initiate_var:
			chrom, position = var.split(':')
			maternal_gt, paternal_gt = get_parental_genotypes(	parents_reader,
										parents_sample_names,
										chrom,
										position)
			if maternal_gt is None and paternal_gt is None:
				continue 
			template_lengths_at_position_dic = get_reads_tlen(bam_reader, chrom, position)
			position_list = []
			initiate_var = False
		line = line.rstrip().split('\t')
		geno = line[3]
		read = line[4].split(':')
		qname = ':'.join(read[1:-10])
		strand = 1 if read[-6] == '+' else 0
		isize = template_lengths_at_position_dic[qname]

		position_list.append([geno, isize, qname, strand])
		# print(position_list)

	elif 'TYPE=' in line:

		line_list = line.rstrip().split('\t')
			
		ref, alt = line_list[3:5]
		one_alt_allele = len(alt.split(',')) == 1

		'''
		if a variant exists in the parental VCF but no reads were found in the cfDNA, no lines will
		start with "haplo_obs", therefore initiate_var will still be set as True. So only write to db
		if initiate_var is set as False.
		'''
		reads_were_found_in_the_cfdna = not initiate_var
		if one_alt_allele and reads_were_found_in_the_cfdna:

			var_type_string = line_list[7].split('TYPE=')[1].split(';')[0]
			var_type = var_type_dic[var_type_string]
			
			format_fields = line_list[8].split(':')

			# print(position_list, file=sys.stderr)
			if args.downsample:
				# F - original fetal fraction, G - new fetal fraction
				F, D = [float(i) for i in args.downsample.split(',')]
				# F = 0.3021516884731206
				# D = 0.13
				fet = which_allele_is_fetal(maternal_gt, paternal_gt, ref, alt)
				position_list = downsample_position(	position_list,
									F,
									D,
									fetal_frequencies,
									maternal_frequencies,
									ref,
									alt,
									fet)
			# print(position_list, file=sys.stderr)
			
			if len(position_list) == 0:
				continue

			for l in position_list:
				genotype = l[0]	
				is_fetal = is_fetal_fragment(	genotype,
								ref,
								alt,
								fetal_allele = get_fetal_allele_type(maternal_gt, paternal_gt))
				for_ff = use_for_fetal_fraction_calculation(	maternal_gt,
										paternal_gt,
										var_type,
										is_fetal)
				l += [is_fetal if is_fetal is not None else 0, for_ff]

			# print(position_list)
			vardb.insertVariant(chrom.replace('chr',''), int(position), position_list)
		
		print(line, end = '')

if args.origin:
	vardb.createQnamesTable()

vardb.lengthDists()

bam_reader.close()

print('finished successfully', file = sys.stderr)