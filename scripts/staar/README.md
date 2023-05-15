---

The pipeline is composed by using [STAAR](https://github.com/xihaoli/STAAR), and the scripts are midified from [STAARpipeline-Tutorial](https://github.com/xihaoli/STAARpipeline-Tutorial).

## Files Description
### Preproperation
- Install xsv [v0.13.0](https://github.com/BurntSushi/xsv)
- Download FAVORdatabase_chrsplit.csv from [STAARpipeline-Tutorial](https://github.com/xihaoli/STAARpipeline-Tutorial/blob/main/FAVORannotator_csv/FAVORdatabase_chrsplit.csv)
- Download FAVORannotator's Full database CSV files
	- [https://favor.genohub.org](https://favor.genohub.org/)
	- Check FAVORannator tab
- Untar FAVORannotator's Full database CSV files

### Input file
- VCF file of a chromosome

### Output file
- gds, agds, and assoc.csv

---

### Execution

1. For each VCF file of a chromosome, generate a GDS [generate_chr_gds.R](generate_chr_gds.R)
2. For each chromosome using the GDS generated from step 1, generate multiple csv files for variants that will be annotated [Varinfo_gds.R](Varinfo_gds.R)
	- Outputs: 
		- A folder `<Out_Path>/chr<CHR>` will be created.
		- `VarInfo_chr$chr_$id.csv` will be gerated in `<Out_Path>/chr<CHR>`
3. Annotate variants of a chromosome [Annotate.R](Annotate.R)
	- Note that the script will search all csv files generated from step 2 and perform annotation
	- Outputs:
		- `Anno_chr$chr_$id.csv`: a CSV file containing annotated variants in `VarInfo_chr$chr_$id.csv`
		- `Anno_chr$chr.csv`: a CSV file contaning all annotated variants in all `Anno_chr$chr_$id.csv`
		- `Anno_chr$chr_STAARpipeline.csv`: a CSV file containing the variants list with annotations required for STAARpipeline of the chromosome.
4. Generate aGDS [gds2agds.R](gds2agds.R) using GDS and `Anno_chr$chr_STAARpipeline.csv`

---

## Dependencies
### Tools
- Biobase
- dplyr
- GENESIS
- GGally
- SeqArray
- SeqVarTools
