
# These gene sets for GSEA analyses were downloaded from http://software.broadinstitute.org/gsea/msigdb/collections.jsp.
# In this page you can find a description for each of them: http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp#C7


# Scripts to perform gtf to bed12 are from UCSC. This is the code used to download them
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod 755 gtfToGenePred genePredToBed
