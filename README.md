# eurl_vtec_wgs_pt-galaxy
 This tool performs various Escherichia coli typing tools and is implemented as a Galaxy (https://galaxyproject.org/) workflow

    Raw data quality check (FASTQC, http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    Trimming (Trimmomatic, DOI:10.1093/bioinformatics/btu170)
    Assembly (SPAdes, DOI:10.1089/cmb.2012.0021; SKESA, DOI:10.1186/s13059-018-1540-z)
    Virulotyping (patho_typing tool from the INNUENDO Project)
    Multi Locus Sequence Typing (MLST 7 loci, https://github.com/tseemann/mlst/)
    Serotyping (MMseqs2, DOI:10.1038/s41467-018-04964-5)
    Shigatoxintyping (blastn, DOI:10.1186/1471-2105-10-421, of a consensus sequence against the shiga toxin subtype database from the Statens Serum Institut SSI and Technical University of Denmark DTU, DOI:10.1128/JCM.00008-15)
    AMR typing (AMRFinderPlus, DOI:10.1128/aac.00483-19)

After installation the BASE_URL parameter in the EURL_VTEC_WGS_PT.py (line 20) will have to be modified to the url of your Galaxy instance in order to correctly visualise the FastQC results.  

The files duk, fastq_pair, stx_subtype_fa.sh, stx_subtype_pe.sh, stx_subtype_se.sh in the scripts folder and the file rematch.py in the scripts/ReMatCh folder should have execution rights.  

In order to make trimmomatic work, you will have to create the following symbolic link
    cd /afs/galaxy/t_d_d/_conda/envs/mulled-v1-9471dd12387e90c11124403c650a667dc2a8c932d610ab6fc4cf2e3f4b40720c/bin
    chmod 755 ../share/trimmomatic-0.39-1/trimmomatic.jar
    ln -s ../share/trimmomatic-0.39-1/trimmomatic.jar trimmomatic.jar
