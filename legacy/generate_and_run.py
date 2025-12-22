import os

# 基础路径
base_path = "../dataset/CAMIMOUSE/output_dir"
file_prefix = "2017.12.29_11.37.26_sample_"
file_suffix_f = "_contigs_2017.12.29_11.37.26_sample_"
file_ext_f = "_contigs_anonymous_gsa.fasta"
file_suffix_sf = "_contigs_2017.12.29_11.37.26_sample_"
file_ext_sf = "_contigs_gsa_mapping_new.tsv"

# 参数部分
software = "ganon ganon2 kraken2 bracken taxor chimera"
database = "complete"
threads = 32
confidence = 0.7

# 构造 -f 参数
f_files = " ".join([f"{base_path}/{file_prefix}{i}{file_suffix_f}{i}{file_ext_f}" for i in range(64)])

# 构造 -sf 参数
sf_files = " ".join([f"{base_path}/{file_prefix}{i}{file_suffix_sf}{i}{file_ext_sf}" for i in range(64)])

# 构造完整指令
command = f"python runClassifier.py -s {software} -d {database} -t {threads} -c {confidence} -f {f_files} -sf {sf_files}"

# 打印生成的指令
print("Generated Command:")
print(command)

# 执行指令
os.system(command)