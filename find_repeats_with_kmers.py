import subprocess
import pandas as pd
import os
from os.path import  exists

# 定义K-mer大小和阈值
kmer_size = 100
Repetitive_quantity = 15
genome_name = "kmar"

# 定义输入FASTA文件和输出Jellyfish文件
genome_fasta_path = "/hpcfs/fhome/yangchh/workdir/self/editSeqDesign/AutoSGRS/input/GCA_001761485.1_ASM176148v1_genomic.fna"
output_path = "/hpcfs/fhome/yangchh/workdir/self/editSeqDesign/AutoSGRS/output"
if not exists(output_path):
    os.makedirs(output_path)  


#拼接前序列
target1_path = os.path.join(output_path, "kmer_100_15.fasta")
#拼接后序列的fasta文件路径
target2_path = os.path.join(output_path, "kmer_100_15_bowtie.fasta")
  
#kmer参数
jellyfish_output = os.path.join(output_path, f"{genome_name}_kmer_counts.jf" )
jellyfish_output_tsv = os.path.join(output_path, f"{genome_name}_kmer_counts.tsv" )

#软件路径
BOWTIE_PATH = "/hpcfs/fhome/yangchh/software/bowtie"
bowtie_workdir = os.path.join(output_path, "bowtie")


def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

def merge_df_by_rule(kmer_350_sam_df):
    df_list = []
    k = 0
    next = 0
    all_list = []
    
    for i,v in kmer_350_sam_df.iterrows():
            
            if next == 0:  
                next = v["ReferenceStart"]     

            start = v["ReferenceStart"]

            if start == next:
                pass
            else:
                print("kljdf")
                df = kmer_350_sam_df[k:i]
                df_list.append(df)
                k = i
                next = start

            next = next + 1

    df = kmer_350_sam_df[k:]
    df_list.append(df)

    for item in df_list:
            str_merge_all = str_merge(list(item["Sequence"]))
            all_list.append(str_merge_all)

    str_merge_df = pd.DataFrame( data = all_list, columns = ["sequence"] ) 
    str_merge_df["id"] = str_merge_df.index
    str_merge_df = str_merge_df.drop_duplicates("sequence")

    return str_merge_df


def str_merge(str_list):

    # 初始化结果字符串为第一个字符串
    merged_result = str_list[0]

    # 从第二个字符串开始，每个字符串只取最后一个字母进行拼接
    for s in str_list[1:]:
        merged_result += s[-1]
    return merged_result

def df_to_fasta(df,fasta_filename):
    with open(fasta_filename, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence_id = row['id']
            sequence = row['sequence']
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')

def parse_sam_file(sam_file_path = "alignment.sam"):

    # 定义用于存储比对结果的列表
    alignment_data = []  

    with open(sam_file_path, "r") as sam_file:
        for line in sam_file:
            if not line.startswith("@"):  # 跳过SAM文件头部
                fields = line.strip().split("\t")
                read_name = fields[0]  
                chain = fields[1]
                reference_name = fields[2]
                reference_start = int(fields[3])  
                sequence = fields[9]
                mismatch = fields[-2]
                matching_number = fields[-1]

                alignment_data.append([read_name, chain, reference_name, reference_start, sequence, mismatch, matching_number])

    # 创建DataFrame
    columns = ["ReadName","chain", "ReferenceName", "ReferenceStart", "Sequence","Mismatch", "MatchingNumber"]
    alignment_df = pd.DataFrame(alignment_data, columns=columns)

    return alignment_df

def bowtie_seq_genome(workdir, genome_path, genome_name, target_path ,option):

    if not exists(workdir):
        os.makedirs(workdir)

    index_path = os.path.join(workdir,genome_name)  

    sam_path = os.path.join( index_path, 'output.sam')  

    if not exists(index_path):
        os.makedirs(index_path)
        # index = os.path.join(index_prefix, 'genome_index')
        # os.chmod(f'{config.BOWTIE_PATH}bowtie-build', 0o755)
        index_prefix = os.path.join(index_path, 'genome_index')
        cmd = f'{BOWTIE_PATH}/bowtie-build {genome_path} {index_prefix}'
        print(cmd)
        os.system(cmd)

    # os.chmod(f'{config.BOWTIE_PATH}bowtie', 0o755)
    index_prefix = os.path.join(index_path, 'genome_index')
    cmd = f'{BOWTIE_PATH}/bowtie -p 2 -v 0  {option} --sam-nohead -k 1000000000 {index_prefix} -f {target_path} -S {sam_path}'
    print(cmd) 
    os.system(cmd)

    #解析  
    sam_df = parse_sam_file(sam_path)  

    #删除文件

    return sam_df


#1. 运行Jellyfish计数  
jellyfish_cmd = f"jellyfish count  -m {kmer_size} -t 16 -s 100M  -o {jellyfish_output} {genome_fasta_path}"  
subprocess.call(jellyfish_cmd, shell=True)
print( jellyfish_cmd )


#2. 按列分别存储信息的计数结果
count_cmd_csv = f"jellyfish dump -c -t {jellyfish_output} > {jellyfish_output_tsv}"
subprocess.call(count_cmd_csv, shell=True)


#3. 读取kmer结果
column_names = ["sequence", "count"]
df = pd.read_csv(jellyfish_output_tsv, sep='\t', names=column_names)
print( "kmer序列的dataframe：" ,df[df["count"]>=Repetitive_quantity] )


#4. 存kmer序列到fasta
kmer_350_df = df[df["count"]>=Repetitive_quantity].reset_index(drop=True) 
kmer_350_df["id"] = kmer_350_df.index.astype(str) + "_" + kmer_350_df["count"].astype(str)
df_to_fasta(kmer_350_df, target1_path)
# 打印第一个 k-mer 序列的长度
if not kmer_350_df.empty:
    print("kmer序列的长度：", len(kmer_350_df.loc[0, "sequence"]))
else:
    print("过滤后的 DataFrame 为空")


#5. 序列比对，拿到kmer序列的坐标
kmer_350_sam_df = bowtie_seq_genome(bowtie_workdir, genome_fasta_path, genome_name, target1_path,option="--norc")
kmer_350_sam_df = kmer_350_sam_df.sort_values(by='ReferenceStart', ascending=True)
kmer_350_sam_df = kmer_350_sam_df.drop_duplicates("Sequence").reset_index(drop=True)


#6. 拼接kmer序列
merge_str_df = merge_df_by_rule(kmer_350_sam_df)


#7. 序列比对拼接后的序列
df_to_fasta(merge_str_df, target2_path)
merge_str_sam_df = bowtie_seq_genome(bowtie_workdir, genome_fasta_path, genome_name, target2_path, option="")
merge_str_sam_df["len"] = merge_str_sam_df.Sequence.apply(lambda x: len(x))         


#8. 对拼接后的序列进行标准化处理
merge_str_sam_df["Sequence"] = merge_str_sam_df.apply(lambda row: revComp(row["Sequence"]) if row["chain"] == "16" else row["Sequence"], axis=1)
merge_str_sam_df["MatchingNumber"] = merge_str_sam_df.MatchingNumber.str.replace("XM:i:", "")
merge_str_sam_df["ReferenceEnd"] = merge_str_sam_df["ReferenceStart"]+merge_str_sam_df["len"]
merge_str_sam_df["Coordinate"] = merge_str_sam_df["ReferenceName"]+":"+merge_str_sam_df["ReferenceStart"].astype(str)+"-"+merge_str_sam_df["ReferenceEnd"].astype(str)
merge_df = merge_str_sam_df[["Sequence", "Coordinate", "len", "MatchingNumber"]].sort_values(by='MatchingNumber', ascending=False)
def sort_coordinates(coords):
    sorted_coords = sorted(coords, key=lambda x: int(x.split(":")[1].split("-")[0]))
    return ",".join(sorted_coords)
merge_df = merge_df.groupby("Sequence")["Coordinate"].agg(lambda coords: sort_coordinates(coords)).reset_index()
merge_df["len"] = merge_df.Sequence.apply(lambda x: len(x))
merge_df["MatchingNumber"] = merge_df.Coordinate.apply(lambda x: len(x.split(",")))

#9. 输出
merge_df["Genome"] = genome_name
merge_df.to_csv(os.path.join(output_path,"result.csv"), index=False)

