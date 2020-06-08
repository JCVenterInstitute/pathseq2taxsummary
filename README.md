# pathseq2taxsummary
Pathseq2taxsummary is a Perl script to convert a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file to a  MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/), which can then be used to make various plots and to compute various statistics using R or other software packages. Another useful feature of this script is that it enables combined analysis of multiple samples while the PathSeq scores.txt file only includes read mapping results from a single sample.

## Examples

Included in this distribution are input and output files that were used in Lang et. al. and Jian et. al. (see Citations).

To run pathseq2taxsummary example datasets, cd into the desired directory (e.g., Lang_dir or Jian_dir) and run the following command(s):

1) Concatenate all score.txt files for each sample in the dataset (appending "{sample_prefix}/scores.txt" at the end of each line using awk) [in bash shell]
 * for i in $(cat sample_prefixes.lst); do awk -v OFS='\t' '{print $0,FILENAME}' $i/scores.txt >> combined_scores.txt; done;

2. Make the .taxsummary file for PathSeq "ambiguous" mapped read counts (see https://gatkforums.broadinstitute.org/gatk/discussion/10913/how-to-run-the-pathseq-pipeline for definitions)

 * ../../pathseq2taxsummery.pl -s combined_scores.txt > combined_scores_taxsummary_amb.txt

3. Make the .taxsummary file for PathSeq "unambiguous" mapped read counts (see https://gatkforums.broadinstitute.org/gatk/discussion/10913/how-to-run-the-pathseq-pipeline for definitions)

 * ../../pathseq2taxsummery.pl -s combined_scores.txt -u > combined_scores_taxsummary_unamb.txt

4. Filter out specific taxa if desired.  In this example, we pulled only "Viruses"
 * head -1 combined_scores_taxsummary_amb.txt > combined_virus_scores_taxsummary_amb.txt 
 * grep "Viruses" combined_scores_taxsummary_amb.txt >> combined_virus_scores_taxsummary_amb.txt
 
 * head -1 combined_scores_taxsummary_unamb.txt > combined_virus_scores_taxsummary_unamb.txt 
 * grep "Viruses" combined_scores_taxsummary_unamb.txt >> combined_virus_scores_taxsummary_unamb.txt

## Usage

**pathseq2taxsummary.pl -s <modified PathSeq scores.txt file> [options]**

	pathseq2taxsummary.pl [options]
	Options:
		-s		: [D]ebug

		-u		: display the last [i]nvocation to the user.

		-h    		: fr[A]gment file [REQUIRED]
				: (<Level><\t><Fragment_id><\t><Left_end><\t><Right_end><\t><Left_margin>)

## Included Files
Two example datasets from Lang et. al. and Jiang et. al.

### pathseq2taxsummary.pl
A PERL script that generates a MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/) from a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file, which can then be used to make various plots and to compute various statistics using R or other software packages.

### examples_dir
A directory containing 3 .tar.bz archives of two directories of example data for generating tax.summary files from Lang et. al. and Jiang et. al.

To unarchive Lang_dir.tar.bz (in the examples_dir):

* cd examples_dir
* tar -xvjf Lang_dir.tar.bz

To unarchive and combined the 2 Jiang et. al. archives (split due to size restrictions) into a single Jiang_dir folder:

* cd examples_dir
* mkdir Jiang_dir
* tar -xvjf Jiang_dir_a.tar.bz -C Jiang_dir --strip-components=1
* tar -xvjf Jiang_dir_b.tar.bz -C Jiang_dir --strip-components=1

## Dependencies

**The following is a list of required Perl modules used by LinearDisplay.pl:**

	Getopt::Std		: Included in Perl 5 distribution
	Cwd			: Included in Perl 5 distribution

## Citation

The following publications used this program to generate linear illustrations of bacteriophage genomes and should be used to site this program:The following publications used this program to convert a slightly modified and concatenated PathSeq (http://software.broadinstitute.org/pathseq/) scores.txt file to a  MOTHUR style tax.summary file (https://mothur.org/wiki/summary.tax/), which can then be used to make various plots and to compute various statistics using R or other software packages.:

Lang S, Demir M, Martin A, Jiang L, Zhang X, Duan Y, Gao B, Wisplinghoff H, Kasper P, Roderburg C, Tacke F, Steffen H-M, Goeser T, Abraldes J G, Tu X M, Loomba R, Pride D, Fouts D E, Schnabl B. Intestinal virome signature associated with non-alcoholic steatohepatitis and fibrosis. Gastroenterology. in press.

Jiang L, Lang S, Duan Y, Zhang X, Gao B, Chopyk J, Schwanemann L K, Ventura-Cots M, Bataller R, Bosques-Padilla F, Verna E C, Abraldes J G, Brown Jr R S, Vargas V, Altamirano J, Caballería J, Shawcross D L, Ho S B, Louvet A, Lucey M R, Mathurin P, Garcia-Tsao G, Kisseleva T, Brenner D A, Tu X M, Stärkel P, Pride D, Fouts D E and Schnabl B. Intestinal virome in patients with alcoholic hepatitis. Hepatology. in press.

## Contact
Derick E. Fouts
dfouts@jcvi.org
