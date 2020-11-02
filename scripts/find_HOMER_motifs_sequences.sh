cd ~/projects/MPRA/WHR/results/motif/enhancers_vs_hg38/

head -39 knownResults.txt | tail -38 | while read line; do
	motif=$(echo $line | cut -f1 | cut -f1 -d'(')
	motif_img=known${i}.logo.svg
	mv knownResults/$motif_img knownResults/${motif}.logo.svg
	motif_file=known${i}.motif
	mv knownResults/$motif_file knownResults/${motif}.motif
	i=$(( $i + 1 ))
done

cd ../../

for i in COUP-TFII Erra THRb Usf2 LXRE EAR2 RARa Atf1; do
	findMotifs.pl motif/fa/EMVars.fa fasta motif/$i -fasta motif/fa/nonsig.fa -find motif/enhancers_vs_hg38/knownResults/$i.motif > motif/EMVars_vs_nonsig/$i.txt
	findMotifs.pl motif/fa/enhancers.fa fasta motif/$i -fasta motif/fa/nonsig.fa -find motif/enhancers_vs_hg38/knownResults/$i.motif > motif/enhancers_vs_nonsig/$i.txt
	findMotifs.pl motif/fa/nonsig.fa fasta motif/$i -fasta motif/fa/nonsig.fa -find motif/enhancers_vs_hg38/knownResults/$i.motif > motif/enhancers_vs_nonsig/${i}_nonsig.txt
done
