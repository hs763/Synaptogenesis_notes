# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

#path to the FastQs
PBS='/data2/hanna/synaptogenesis'

PATH="/data2/hanna/synaptogenesis/ParseBiosciences-Pipeline.1.0.3p:$PATH"

cd $PBS/newvolume/genomes/
#cat Homo_sapiens.GRCh38.108.gtf EmGFP.gtf > Homo_sapiens.GRCh38.108.EmGFP.gtf
#cat Homo_sapiens.GRCh38.dna.primary_assembly.fa EmGFP.fa > Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa

#cat Mus_musculus.GRCm39.108.gtf EmGFP.gtf > Mus_musculus.GRCm38.108.EmGFP.gtf
#cat Mus_musculus.GRCm39.dna.primary_assembly.fa EmGFP.fa > Mus_musculus.GRCm38.dna.primary_assembly.EmGFP.fa

#split-pipe \
#--mode mkref \
#--genome_name hg38 \
#--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa \
#--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.EmGFP.gtf \
#--output_dir $PBS/newvolume/genomes/hg38

split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.gtf \
--output_dir $PBS/newvolume/genomes/hg38_noGFP

#cat Homo_sapiens.GRCh38.108.gtf hEmGFP.gtf > Homo_sapiens.GRCh38.108.hEmGFP.gtf
#cat Homo_sapiens.GRCh38.dna.primary_assembly.fa hEmGFP.fa > Homo_sapiens.GRCh38.dna.primary_assembly.hEmGFP.fa

#cat Mus_musculus.GRCm39.108.gtf mEmGFP.gtf > Mus_musculus.GRCm38.108.mEmGFP.gtf
#cat Mus_musculus.GRCm39.dna.primary_assembly.fa mEmGFP.fa > Mus_musculus.GRCm38.dna.primary_assembly.mEmGFP.fa

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 mm10 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa $PBS/newvolume/genomes/Mus_musculus.GRCm38.dna.primary_assembly.fa \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.gtf $PBS/newvolume/genomes/Mus_musculus.GRCm38.108.gtf \
--output_dir $PBS/newvolume/genomes/hg38_mm10_noGFP

cd $PBS/newvolume/expdata/

cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_1.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_1.fq.gz > SLX-22602.tmp.r_1.fq.gz
cat SLX-22602.DNAA007.HGMLNDMXY.s_1.r_2.fq.gz SLX-22602.DNAA007.HGMLNDMXY.s_2.r_2.fq.gz > SLX-22602.tmp.r_2.fq.gz

cat SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz > SLX-22602.lostreads.tmp.r_1.fq.gz
cat SLX-22602.HGMLNDMXY.s_2.r_1.lostreads.fq.gz SLX-22602.HGMLNDMXY.s_2.r_2.lostreads.fq.gz > SLX-22602.lostreads.tmp.r_2.fq.gz

cat SLX-22602.tmp.r_1.fq.gz SLX-22602.lostreads.tmp.r_1.fq.gz > SLX-22602.r_1.fq.gz
cat SLX-22602.tmp.r_2.fq.gz SLX-22602.lostreads.tmp.r_2.fq.gz > SLX-22602.r_2.fq.gz

#not run currently
nohup ./demultiplexer.rhel/demuxFQ \
    -c -d -e -i -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22602.lostreads.r_1.fq.gz \
    -s SLX-22602.demultiplexsummary.r1.txt \
    SLX-22602.r_1.index.txt \
    SLX-22602.r_1.fq.gz &
    
nohup ./demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22602.lostreads.r_2.fq.gz \
    -s SLX-22602.demultiplexsummary.r2.txt \
    SLX-22602.r_2.index.txt \
    SLX-22602.r_2.fq.gz &
    
rm SLX-22602.r_1.fq.gz SLX-22602.r_2.fq.gz



# Pipeline running
#single cell
#for mixed specie genome 
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/ACTTGA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/AGTCAA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/AGTTCC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/ATGTCA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/CAGATC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/CTTGTA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/GATCAG

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38_mm10_noGFP/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell/TAGCTT

split-pipe \
    --mode comb \
    --sublibraries $PBS/newvolume/analysis/sCell/ACTTGA $PBS/newvolume/analysis/sCell/AGTCAA $PBS/newvolume/analysis/sCell/AGTTCC $PBS/newvolume/analysis/sCell/ATGTCA $PBS/newvolume/analysis/sCell/CAGATC $PBS/newvolume/analysis/sCell/CTTGTA $PBS/newvolume/analysis/sCell/GATCAG  $PBS/newvolume/analysis/sCell/TAGCTT \
    --output_dir $PBS/newvolume/analysis/sCell/combined


#human reference genome
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/ACTTGA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/AGTCAA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/AGTTCC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/ATGTCA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/CAGATC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/CTTGTA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/GATCAG

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/hg38/TAGCTT &

split-pipe \
    --mode comb \
    --sublibraries $PBS/newvolume/analysis/hg38/ACTTGA $PBS/newvolume/analysis/hg38/AGTCAA $PBS/newvolume/analysis/hg38/AGTTCC $PBS/newvolume/analysis/hg38/ATGTCA $PBS/newvolume/analysis/hg38/CAGATC $PBS/newvolume/analysis/hg38/CTTGTA $PBS/newvolume/analysis/hg38/GATCAG  $PBS/newvolume/analysis/hg38/TAGCTT \
    --output_dir $PBS/newvolume/analysis/hg38/combined_hg38
    
    
    
    
#mouse reference genome
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/ACTTGA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/AGTCAA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/AGTTCC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/ATGTCA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/CAGATC

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/CTTGTA

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/GATCAG

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/mm10/ \
--fq1 $/data2/ivanir/Feline2023/ParseBS/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $For human brian organoids: 
The EB medium has been changed on day 3, where I removed as much of the old medium without loosing the organoid, washed once with 150 ul of EB medium (without rock inhibitor), and replaced with 150 ul of EB medium. On day 5, I replaced the EB medium with 150 ul neural induction (NI) medium, with a washing step as before. On day 7, the organoids were embedded in matrigel droplet. I cut the parafilm into squares, sprayed with ethanol, indented in 16 places and left to try. I transferred the organoids form wells to the indentations using cut 200 ul pipette set to 150 ul, and removed excess media. Add a drop of thawed matrigel aliquot ont each organoid and use 20 ul pipet tip to surround the organoid with matrigel from all sided. I incubated the matrigel coated organoids at 37°C, 5% CO2 for 20 min. I then transferred the organoids into wells of a 6 well plate filled with 3 ml of expansion medium, with max. 16 organoids per well, and incubated at 37°C, 5% CO2. On day 6, I changed the medium into Maturation medium. On day 10, medium was changed to IDM -A. On day 14, I removed the matrigel initially using needles under a microscope to cut mechanically remove most of it. Then, I transferred max. 6 organoids at the time into a sterile eppendorf tube using cut 1 ml pipet tip and removed all the excess media. I incubated them in 500 ml of rCell REcovery Solution (CLS354253, Corning)  solution at 4°C ,for 20 minutes. I removed the recovery solution, washed three times with PBS, and suspended them in 1 ml of warm IDM -A + matrigel. I incubated the organoids for 2 h at 37°C, 5% without orbital shaker, making sure organoids are not touching to prevent fusion. I den transferred the petri dish with organoids to orbital shaker at 37°C, 5% CO2. The organoids were fed 2x a week with IDM-A for the frost week and with IDM+A + matrigel for the second week. 
/newvolume/expdata/correctFastq/SLX-22602.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/mm10/TAGCTT

split-pipe \
    --mode comb \
    --sublibraries $PBS/newvolume/analysis/mm10/ACTTGA $PBS/newvolume/analysis/mm10/AGTCAA $PBS/newvolume/analysis/mm10/AGTTCC $PBS/newvolume/analysis/mm10/ATGTCA $PBS/newvolume/analysis/mm10/CAGATC $PBS/newvolume/analysis/mm10/CTTGTA $PBS/newvolume/analysis/mm10/GATCAG  $PBS/newvolume/analysis/mm10/TAGCTT \
    --output_dir $PBS/newvolume/analysis/sCell/combined_mm10
