# this is a protocol that we used in
# "Common Structural Patterns in the Maxicircle Divergent Region of Trypanosomatidae"
# in MDPI Pathogens (https://doi.org/10.3390/pathogens9020100)
# to assemble full-length maxicircle sequences


# exports for our system
#MUMMER_DIR=/home/jalgard/Tools/MUMmer3.23
#TOOL_DIR=/home/jalgard/Tools
#CANU_DIR=/home/jalgard/Tools/Canu/bin

#download SRR
SRR={srr}

# for further usage with Canu it will go (despite the reads are actually from PacBio)
fastq-dump -I $SRR

# convert to fasta
cat SRR*.fastq | awk 'BEGIN{rc=0;}{if(NR % 4 == 1) {rc=rc+1; print ">" rc;} if(NR % 4 == 2) {print $0}}' > $SRR.fasta

#make blastdb from all pacbio reads
makeblastdb -in $SRR.fasta -out $SRR.db -dbtype nucl

#find reads with 12S and ND5 traces
blastn -query ../12S.nuc.fasta -db $SRR.db -evalue 1e-10 -outfmt 6 -num_threads 5 -word_size 7 > $SRR.12S_blast_out
blastn -query ../ND5.nuc.fasta -db $SRR.db -evalue 1e-10 -outfmt 6 -num_threads 5 -word_size 7 > $SRR.ND5_blast_out

# join reads
cat $SRR.12S_blast_out $SRR.ND5_blast_out | cut -f2 | sort | uniq > $SRR.maxicircle.reads
# pick (this stage uses custom python script available in my repository 'totoro')
python $TOOL_DIR/fasta-kit.py --in $SRR.fasta --list $SRR.maxicircle.reads --action keep --out $SRR.maxi.selected.fasta

# assemble with Canu (Zk - genome size, we used Z=50 for most assemblies)
$CANU_DIR/canu -pacbio-raw $SRR.maxi.selected.fasta -d $SRR.canu -p $SRR genomeSize=Zk obtovlThreads=20 merylThreads=20 corThreads=20 ovlThreads=20 corMaxEvidenceErate=0.15 correctedErrorRate=0.105

# check for circle
$MUMMER_DIR/nucmer -maxmatch -nosimplify $SRR.canu/$SRR.contigs.fasta $SRR.canu/$SRR.contigs.fasta
$MUMMER_DIR/show-coords -lrcTH out.delta > $SRR.assm.delta

# use mreps to find tandem repeats
$TOOL_DIR/mreps/mreps -res 5 -minsize 10 -minperiod 4 -fasta $SRR.canu/$SRR.contigs.fasta > $SRR.mreps
# find inverted repeats
$TOOL_DIR/emboss/bin/einverted -sequence $SRR.contigs.fasta -outfile $SRR.einverted -outseq $SRR.invseq -match 3 -mismatch -4 -gap 12 -threshold 50
