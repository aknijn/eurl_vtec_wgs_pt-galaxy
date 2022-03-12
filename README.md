# eurl_vtec_wgs_pt-galaxy
 This tool performs various Escherichia coli typing tools and is implemented as a Galaxy (https://galaxyproject.org/) workflow

    Raw data quality check (FASTQC, http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    Trimming (Trimmomatic, DOI:10.1093/bioinformatics/btu170)
    Assembly (SPAdes, DOI:10.1089/cmb.2012.0021; SKESA, DOI:10.1186/s13059-018-1540-z)
    Virulotyping (patho_typing tool from the INNUENDO Project)
    Multi Locus Sequence Typing (MLST 7 loci, Seemann T, https://github.com/tseemann/mlst/)
    Shigatoxintyping (blastn, DOI:10.1186/1471-2105-10-421, of a consensus sequence against the shiga toxin subtype database from the Statens Serum Institut SSI and Technical University of Denmark DTU, DOI:10.1128/JCM.00008-15)
    AMR typing (Abricate, Seemann T, Github https://github.com/tseemann/abricate; ResFinder database DOI:10.1093/jac/dks261)

After installation the BASE_URL parameter in the EURL_VTEC_WGS_PT.py file (line 20) will have to be modified to the url of your Galaxy instance in order to correctly visualise the FastQC results.  
