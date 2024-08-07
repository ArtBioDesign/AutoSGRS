  
# AutoSGRS
## Project Introduction  
**AutoSGRS**is a tool that automatically searches for repeated sequences of any frequency in any genome

![AutoSGRS]() 
## Installation
### python packages
We suggest using Python 3.8 for AutoSGRS.
```shell
pip install -r requirements.txt

```

### bowtie
```shell
wget https://bowtie-x64-linux.tar.gz -O ~/bowtie-x64-linux.tar.gz
tar -zxvf bowtie-x64-linux.tar.gz
export PATH=~/bowtie/bin:$PATH




## Usage & Example

**Input:**
- **Step 1:** Upload destination genome file(fasta)
- **Step 2:** provide the necessary configuration information.
    - Example configuration:
      ```json
      {
        kmer_size = 100
        Repetitive_quantity = 15
        genome_name = "Yarrowia_lipolytica"
      }
      ```   
      

**Execute:**

```shell
python find_repeats_with_kmers.py
```
**Output:**
- `xxx.csv` 
These files will be generated in the `XXX/AutoSGRS/output` directory.

