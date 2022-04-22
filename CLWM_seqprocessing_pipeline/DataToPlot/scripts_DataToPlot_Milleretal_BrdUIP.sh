###CORRELATION###
##HU samples
#plot correlations
multiBigwigSummary bins -b DataToPlot/bigwig/SPT6-2-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/SPT6-3-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt650-1-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt650-4-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt6YW-2-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt6YW-1-HU-IP_spikeinNODUP.bw -out DataToPlot/spikeinNODUP_forMilleretal_multibigwigsummary.npz --labels WT1 WT2 s501 s502 YW1 YW2

plotCorrelation --corData DataToPlot/spikeinNODUP_forMilleretal_multibigwigsummary.npz --corMethod pearson --whatToPlot heatmap -o DataToPlot/spikeinNODUP_HUMilleretal_pearsonheatmap.png --plotNumbers

###computematrix###
##HU
#generate averaged bigwig
bigwigCompare -b1 DataToPlot/bigwig/SPT6-2-HU-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/SPT6-3-HU-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/SPT6-23avg-HU-IP_spikeinNODUP.bw -of bigwig

bigwigCompare -b1 DataToPlot/bigwig/spt6YW-1-HU-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/spt6YW-2-HU-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/spt6YW-12avg-HU-IP_spikeinNODUP.bw -of bigwig

bigwigCompare -b1 DataToPlot/bigwig/spt650-1-HU-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/spt650-4-HU-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/spt650-14avg-HU-IP_spikeinNODUP.bw -of bigwig

#calculate compute matrix
computeMatrix reference-point -S DataToPlot/bigwig/SPT6-23avg-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt6YW-12avg-HU-IP_spikeinNODUP.bw DataToPlot/bigwig/spt650-14avg-HU-IP_spikeinNODUP.bw -R /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.bed -a 8000 -b 8000 --referencePoint center --smartLabels -o DataToPlot/computemat/computemat_HUaverage_WT23YW12s5014_Nied.mat --outFileNameMatrix DataToPlot/computemat/computemat_HUaverage_WT23YW12s5014_Nied.tab

#plotheatmap
plotHeatmap -m DataToPlot/computemat/computemat_HUaverage_WT23YW12s5014_Nied.mat --colorMap RdYlBu_r --dpi 350 --sortUsingSamples 1 --refPointLabel=ARS --perGroup --samplesLabel WT spt6YW spt650 -out DataToPlot/computemat/plotheatmap_HUaverage_NiedARS8kB.png


##30min
#generate average bigwig
bigwigCompare -b1 DataToPlot/bigwig/SPT6-2-30-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/SPT6-3-30-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/SPT6-23avg-30-IP_spikeinNODUP.bw -of bigwig

bigwigCompare -b1 DataToPlot/bigwig/spt6YW-2-30-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/spt6YW-3-30-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/spt6YW-23avg-30-IP_spikeinNODUP.bw -of bigwig

bigwigCompare -b1 DataToPlot/bigwig/spt650-1-30-IP_spikeinNODUP.bw -b2 DataToPlot/bigwig/spt650-4-30-IP_spikeinNODUP.bw --operation mean -o DataToPlot/bigwig/spt650-14avg-30-IP_spikeinNODUP.bw -of bigwig

#calculate the compute matrix 
computeMatrix reference-point -S DataToPlot/bigwig/SPT6-23avg-30-IP_spikeinNODUP.bw DataToPlot/bigwig/spt6YW-23avg-30-IP_spikeinNODUP.bw DataToPlot/bigwig/spt650-14avg-30-IP_spikeinNODUP.bw -R /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.bed -a 8000 -b 8000 --referencePoint center --smartLabels -o DataToPlot/computemat/computemat_30average_WT23YW23s5014_Nied.mat --outFileNameMatrix DataToPlot/computemat/computemat_30average_WT23YW23s5014_Nied.tab

#plotHeatmap 
plotHeatmap -m DataToPlot/computemat/computemat_30average_WT23YW23s5014_Nied.mat --colorMap RdYlBu_r --dpi 350 --sortUsingSamples 1 --refPointLabel=ARS --perGroup --samplesLabel WT spt6YW spt650 -out DataToPlot/computemat/plotheatmap_30average_NiedARS8kB.png

###PEAK CALLING###
##HU samples
#run macs2 on individual samples
macs2 callpeak -t fastq/aligned/SPT6-2-HU-IP_Scer.bam -c fastq/aligned/SPT6-2-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-3-HU_q.01_ext300 --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 
macs2 callpeak -t fastq/aligned/SPT6-3-HU-IP_Scer.bam -c fastq/aligned/SPT6-3-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-3-HU_q.01_ext300 --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt6YW-1-HU-IP_Scer.bam -c fastq/aligned/spt6YW-1-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt6YW-1-HU_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt6YW-2-HU-IP_Scer.bam -c fastq/aligned/spt6YW-2-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt6YW-2-HU_q.01_ext300 --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt650-1-HU-IP_Scer.bam -c fastq/aligned/spt650-1-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-1-HU_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt650-2-HU-IP_Scer.bam -c fastq/aligned/spt650-2-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-2-HU_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt650-4-HU-IP_Scer.bam -c fastq/aligned/spt650-4-HU-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-4-HU_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

#sort peak calling
sort -k8,8nr DataToPlot/macs2/Milleretal/WT-2-HU_q.01_ext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/WT-2-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak

sort -k8,8nr DataToPlot/macs2/Milleretal/WT-3-HU_q.01_ext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/WT-3-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak
sort -k8,8nr DataToPlot/macs2/Milleretal/spt6YW-1-HU_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt6YW-1-HU_q.01_nomodelext300_peaks.sorted.narrowPeak

sort -k8,8nr DataToPlot/macs2/Milleretal/spt6YW-2-HU_q.01_ext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt6YW-2-HU_q.01_ext300_peaks.sorted.narrowPeak
sort -k8,8nr DataToPlot/macs2/Milleretal/spt650-1-HU_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt650-1-HU_q.01_nomodelext300_peaks.sorted.narrowPeak

sort -k8,8nr DataToPlot/macs2/Milleretal/spt650-2-HU_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt650-2-HU_q.01_nomodelext300_peaks.sorted.narrowPeak

#run idr analysis
idr --samples DataToPlot/macs2/Milleretal/WT-2-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak DataToPlot/macs2/Milleretal/WT-3-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file WT-23-idr --plot --log-output-file WT.idr.log

idr --samples DataToPlot/macs2/Milleretal/WT-2-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak DataToPlot/macs2/Milleretal/WT-3-HU_q.01_ext300_peaks.narrowPeak.sorted.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file DataToPlot/macs2/Milleretal/WT-23-idrsignal --plot --log-output-file DataToPlot/macs2/Milleretal/WT.idrsignal.log

idr --samples DataToPlot/macs2/Milleretal/spt6YW-1-HU_q.01_nomodelext300_peaks.sorted.narrowPeak DataToPlot/macs2/Milleretal/spt6YW-2-HU_q.01_ext300_peaks.sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file DataToPlot/macs2/Milleretal/YW-12-idr --plot --log-output-file DataToPlot/macs2/Milleretal/YW.idr.log

idr --samples DataToPlot/macs2/Milleretal/spt6YW-3-HU_q.01_nomodelext300_peaks.sorted.narrowPeak DataToPlot/macs2/Milleretal/spt6YW-2-HU_q.01_ext300_peaks.sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file DataToPlot/macs2/Milleretal/YW-32-idr --plot --log-output-file DataToPlot/macs2/Milleretal/YW23.idr.log

idr --samples DataToPlot/macs2/Milleretal/spt650-1-HU_q.01_nomodelext300_peaks.sorted.narrowPeak DataToPlot/macs2/Milleretal/spt650-2-HU_q.01_nomodelext300_peaks.sorted.narrowPeak --input-file-type narrowPeak --rank p.value --output-file DataToPlot/macs2/Milleretal/s50-12-idr --plot --log-output-file DataToPlot/macs2/Milleretal/s50.idr.log

idr --samples DataToPlot/macs2/Milleretal/spt650-1-HU_q.01_nomodelext300_peaks.sorted.narrowPeak DataToPlot/macs2/Milleretal/spt650-2-HU_q.01_nomodelext300_peaks.sorted.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file DataToPlot/macs2/Milleretal/s50-12-idrsignal --plot --log-output-file DataToPlot/macs2/Milleretal/s50.idrsignal.log

#filter idr > 1kB
awk '{if($3-$2 >= 1000) print}' DataToPlot/macs2/Milleretal/WT-23-idr.bed > DataToPlot/macs2/Milleretal/WT-23-idr.1kB.bed
awk '{if($3-$2 >= 1000) print}' DataToPlot/macs2/Milleretal/YW-12-idr.bed > DataToPlot/macs2/Milleretal/YW-12-idr.1kB.bed
awk '{if($3-$2 >= 1000) print}' DataToPlot/macs2/Milleretal/s50-12-idr.bed > DataToPlot/macs2/Milleretal/s50-12-idr.1kB.bed
#intersect idr results with NiedARS file /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed
bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/WT-23-idr.1kB.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.bed 

bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/YW-12-idr.1kB.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.bed 

bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/s50-12-idr.1kB.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.bed 

#remove duplicates in ARS overlap file

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.bed

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.sorted.bed

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.sorted.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.merged.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.sorted.merged.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.sorted.merged.bed

#intersect WT vs spt6 mutants

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.sorted.merged.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvsYW12idr.sorted.merged.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.sorted.merged.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvss50124idr.sorted.merged.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW12idr.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.merged.bed -v > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvsYW12idr.YWONLY.sorted.merged.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s50124idr.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idr.sorted.merged.bed -v > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvss50124idr.s50ONLY.sorted.merged.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvsYW12idr.YWONLY.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvss50124idr.s50ONLY.sorted.merged.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YWONLYvss50ONLY.sorted.merged.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvsYW12idr.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT23idrvss50124idr.sorted.merged.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YWoverlapvss50overlap.sorted.merged.bed


##30minute samples
##macs2 on individual samples

macs2 callpeak -t fastq/aligned/SPT6-1-30-IP_Scer.bam -c fastq/aligned/SPT6-1-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-1-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue

macs2 callpeak -t fastq/aligned/SPT6-2-30-IP_Scer.bam -c fastq/aligned/SPT6-2-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-2-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue

macs2 callpeak -t fastq/aligned/SPT6-3-30-IP_Scer.bam -c fastq/aligned/SPT6-3-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-3-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt6YW-2-30-IP_Scer.bam -c fastq/aligned/spt6YW-2-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt6YW-2-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt6YW-3-30-IP_Scer.bam -c fastq/aligned/spt6YW-3-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt6YW-3-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt650-1-30-IP_Scer.bam -c fastq/aligned/spt650-1-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-1-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt6YW-1-30-IP_Scer.bam fastq/aligned/spt6YW-3-30-IP_Scer.bam -c fastq/aligned/spt6YW-1-30-input_Scer.bam fastq/aligned/spt6YW-3-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt6YW-13-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak -t fastq/aligned/spt650-4-30-IP_Scer.bam -c fastq/aligned/spt650-4-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-4-30_q.01_nomodelext300 --nomodel --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak --broad -t fastq/aligned/SPT6-1-30-IP_Scer.bam -c fastq/aligned/SPT6-1-30-input_Scer.bam -f BAM -g 11.5e6 --bdg -q 0.01 -n WT-1-30_q.01_broadext300 --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

macs2 callpeak --broad -t fastq/aligned/spt650-1-30-IP_Scer.bam -c fastq/aligned/spt650-1-30-input_Scer.bam --nomodel -f BAM -g 11.5e6 --bdg -q 0.01 -n spt650-1-30_q.01_nomodelbroadext300 --extsize 300 --outdir DataToPlot/macs2/playingwithqvalue 

#merge peaks by 5kB intervals
sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/WT-1-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/WT-1-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/WT-1-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/WT-1-30_q.01_nomodelext300_peaks.5kBmerge.bed

sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/WT-2-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/WT-2-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/WT-2-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/WT-2-30_q.01_nomodelext300_peaks.5kBmerge.bed

sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/spt6YW-13-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/YW-13-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/YW-13-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/YW-13-30_q.01_nomodelext300_peaks.5kBmerge.bed

sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/spt6YW-2-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/YW-2-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/YW-2-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/YW-2-30_q.01_nomodelext300_peaks.5kBmerge.bed

sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/spt650-1-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt650-1-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/spt650-1-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/spt650-1-30_q.01_nomodelext300_peaks.5kBmerge.narrowPeak

sort -k 1,1 -k2,2n DataToPlot/macs2/playingwithqvalue/spt650-4-30_q.01_nomodelext300_peaks.narrowPeak > DataToPlot/macs2/Milleretal/spt650-4-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak

bedtools merge -i DataToPlot/macs2/Milleretal/spt650-4-30_q.01_nomodelext300_peaks.narrowPeak.sortedcoord.narrowPeak -d 5000 > DataToPlot/macs2/Milleretal/spt650-4-30_q.01_nomodelext300_peaks.5kBmerge.narrowPeak
#intersect replicates
bedtools intersect -a DataToPlot/macs2/Milleretal/WT-1-30_q.01_nomodelext300_peaks.5kBmerge.bed -b DataToPlot/macs2/Milleretal/WT-2-30_q.01_nomodelext300_peaks.5kBmerge.bed > DataToPlot/macs2/Milleretal/WT-12-30.5kBmergeINTERSECT.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/YW-13-30_q.01_nomodelext300_peaks.5kBmerge.bed -b DataToPlot/macs2/Milleretal/YW-2-30_q.01_nomodelext300_peaks.5kBmerge.bed > DataToPlot/macs2/Milleretal/YW-13.2-30.5kBmergeINTERSECT.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/spt650-1-30_q.01_nomodelext300_peaks.5kBmerge.narrowPeak -b DataToPlot/macs2/Milleretal/spt650-4-30_q.01_nomodelext300_peaks.5kBmerge.narrowPeak > DataToPlot/macs2/Milleretal/s50-14-30.5kBmergeINTERSECT.bed

#intersect merge results with NiedARS file /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed
bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/WT-12-30.5kBmergeINTERSECT.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.bed 

bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/YW-13.2-30.5kBmergeINTERSECT.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.bed 

bedtools intersect -a /n/groups/winston/clw20/0_ref/genomes/BEDFILES/Nieduszynski_ARS_Scerheading.slop1kB.bed -b DataToPlot/macs2/Milleretal/s50-14-30.5kBmergeINTERSECT.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.bed 

#remove duplicates in ARS overlap file

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.bed

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.sorted.bed

sort -k1,1 -k2,2n DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.sorted.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.merged.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.sorted.merged.bed

bedtools merge -i DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.sorted.bed > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.sorted.merged.bed

#compare WT vs spt6 mutants

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.merged.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.YW13.30.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_YW13.30.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.merged.bed -v > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.YW13.YWonly.30.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.merged.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.s5014.30.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_s5014.30.sorted.merged.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.30.sorted.merged.bed -v > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.s5014.s50ONLY.30.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.YW13.YWonly.30.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.s5014.s50ONLY.30.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.spt6.spt6overlap.30.bed

bedtools intersect -a DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.YW13.30.bed -b DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.s5014.30.bed -wa > DataToPlot/macs2/Milleretal/ARSintersect_Nied1kBslop_WT12.vs.spt6.ALLoverlap.30.bed
