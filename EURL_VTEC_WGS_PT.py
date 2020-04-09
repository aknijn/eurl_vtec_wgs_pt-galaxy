#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import sys
import os
import HTML
import datetime
import fileinput

BASE_URL = 'https://aries.iss.it'
TOOL_DIR = os.path.dirname(os.path.abspath(__file__))

def insertFile(filename, report):
    with open(filename) as html_in:
        for line in html_in:
            report.write(line)

def insertFileAsTable(filename, report, hasheader=False, tabclass="table table-rep"):
    with open(filename) as table_in:
        table_data = [[str(col) for col in row.split('\t')] for row in table_in]
    insertTable(table_data, report, hasheader, tabclass)

def insertTable(table_data, report, hasheader=False, tabclass="table table-rep"):
    if hasheader:
        htmlcode = HTML.table(table_data[1:], attribs={'class':tabclass}, header_row=table_data[0])
    else:
        htmlcode = HTML.table(table_data, attribs={'class':tabclass})
    report.write(htmlcode)

def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--serotyping', dest='serotyping', help='perform serotyping', action='store_true')
    parser.add_argument('--virulotyping', dest='virulotyping', help='perform virulotyping', action='store_true')
    parser.add_argument('--shigatoxintyping', dest='shigatoxintyping', help='perform shigatoxintyping', action='store_true')
    parser.add_argument('--amrtyping', dest='amrtyping', help='perform amrtyping', action='store_true')
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input1_ext', dest='input1_ext', help='extension of forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input1_name', dest='input1_name', help='name of forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('-2', '--input2', dest='input2', help='reverse reads file in Sanger FASTQ format')
    parser.add_argument('--input2_ext', dest='input2_ext', help='extension of reverse reads file in Sanger FASTQ format')
    parser.add_argument('--input2_name', dest='input2_name', help='name of reverse reads file in Sanger FASTQ format')
    parser.add_argument('--html1', dest='html1', help='html FASTQC file')
    parser.add_argument('--html1_id', dest='html1_id', help='html FASTQC file id')
    parser.add_argument('--html1_path', dest='html1_path', help='html FASTQC file path')
    parser.add_argument('--text1', dest='text1', help='text FASTQC file')
    parser.add_argument('--html2', dest='html2', help='html FASTQC file')
    parser.add_argument('--html2_id', dest='html2_id', help='html FASTQC file id')
    parser.add_argument('--html2_path', dest='html2_path', help='html FASTQC file path')
    parser.add_argument('--text2', dest='text2', help='text FASTQC file')
    parser.add_argument('--log', dest='logfile', help='log file')
    parser.add_argument('--virulotyper', dest='virulotyper', help='Virulotyping Mapping reads')
    parser.add_argument('--virulotyper_id', dest='virulotyper_id', help='Virulotyping Mapping reads id')
    parser.add_argument('--stx', dest='stx', help='Shiga toxin')
    parser.add_argument('--mlstsevenloci', dest='mlstsevenloci', help='Multi Locus Alleles table')
    parser.add_argument('--amr', dest='amr', help='SPAdes log')
    parser.add_argument('--amr_id', dest='amr_id', help='AMR file id')
    parser.add_argument('--antigen_O', dest='antigen_O', help='Antigen for O')
    parser.add_argument('--antigen_H', dest='antigen_H', help='Antigen for H')
    parser.add_argument('--output', dest='output', help='output report html file')
    args = parser.parse_args()

    log = open(args.logfile, 'w')
    log.write("EURL VTEC WGS PT v3.0\n\nTool versions\n=============\n")
    if args.input2:
        # FASTQC
        os.system("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input1 + " -d " + args.html1_path + " -o " + args.html1 + " -t " + args.text1 + " -f " + args.input1_ext + " -j " + args.input1_name + " -e " + "fastqc")
        os.system("rm -r " + args.html1_path)
        os.system("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input2 + " -d " + args.html2_path + " -o " + args.html2 + " -t " + args.text2 + " -f " + args.input2_ext + " -j " + args.input2_name + " -e " + "fastqc")
        os.system("rm -r " + args.html2_path)
        log.write(os.popen("fastqc -v").read())
        # TRIMMING
        os.system("python " + TOOL_DIR + "/scripts/fastq_positional_quality_trimming.py -1 " + args.input1 + " --maxlt 300 --lt 17 --rt 0 --minqt 25 --avgqt 27.0 --minlf -1 --trimmed1 input_t1.fq --log trimming_logfile -2 " + args.input2 + " --trimmed2 input_t2.fq --trimmedunpaired trimmedunpaired")
        log.write("\nfastq_positional_quality_trimming v1.0\n")
        log.write("parameters: maxlt=300, lt=17, rt=0, minqt=25, avgqt=27.0, minlf=-1\n")
        if args.shigatoxintyping or args.serotyping:
            # ASSEMBLY
            os.system("perl " + TOOL_DIR + "/scripts/spades.pl contigs.fa spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t \${GALAXY_SLOTS:-16} --pe1-ff --pe1-1 fastq:input_t1.fq --pe1-2 fastq:input_t2.fq")
            log.write(os.popen("spades.py -v").read())
            log.write("parameters: --isolate, pe1-ff, pe1-1, pe1-2\n\n")
        if args.virulotyping:
            # VIRULOTYPER
            os.system("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f " + args.input1 + " " + args.input2 + " -o output_dir -j 1 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'")
            os.system("cat pathotyper_rep_tot_tab > " + args.virulotyper)
            log.write("\n\nViruloTyper\n===========\npatho_typing v1.0\n")
            log.write("parameters: minGeneCoverage=90, minGeneIdentity=90, minGeneDepth=15\n\n")
            log.write(os.popen("cat " +  TOOL_DIR + "/data/ViruloTyping_db.txt").read())
    else:
        # FASTQC
        os.system("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input1 + " -d " + args.html1_path + " -o " + args.html1 + " -t " + args.text1 + " -f " + args.input1_ext + " -j " + args.input1_name + " -e " + "fastqc")
        os.system("rm -r " + args.html1_path)
        log.write(os.popen("fastqc -v").read())
        # TRIMMING
        os.system("python " + TOOL_DIR + "/scripts/fastq_positional_quality_trimming.py -1 " + args.input1 + " --maxlt 360 --lt 10 --rt 0 --minqt 25 --avgqt 27.0 --minlf 50 --trimmed1 input_t1.fq --log trimming_logfile")
        log.write("\nfastq_positional_quality_trimming v1.0\n")
        log.write("parameters: maxlt=360, lt=10, rt=0, minqt=25, avgqt=27.0, minlf=50\n")
        if args.shigatoxintyping or args.serotyping:
            # ASSEMBLY
            os.system("perl " + TOOL_DIR + "/scripts/spades.pl contigs.fa spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} --iontorrent -s fastq:input_t1.fq")
            log.write(os.popen("spades.py -v").read())
            log.write("parameters: --isolate, --iontorrent\n\n")
        if args.virulotyping:
            # VIRULOTYPER
            os.system("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f " + args.input1 + " -o output_dir -j 1 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'")
            os.system("cat pathotyper_rep_tot_tab > " + args.virulotyper)
            log.write("\n\nViruloTyper\n===========\npatho_typing v1.0\n")
            log.write("parameters: minGeneCoverage=90, minGeneIdentity=90, minGeneDepth=15\n\n")
            log.write(os.popen("cat " +  TOOL_DIR + "/data/ViruloTyping_db.txt").read())
    # SEQUENCETYPER
    if args.input2:
        os.system("mentalist call --output_votes -o 'mentalist_out' --db '" + TOOL_DIR + "/data/escherichia_coli_pubmlst_k31_2018-10-09/escherichia_coli_pubmlst_k31_m023_2018-10-09.jld' -1 input_t1.fq -2 input_t2.fq")
    else:
        os.system("mentalist call --output_votes -o 'mentalist_out' --db '" + TOOL_DIR + "/data/escherichia_coli_pubmlst_k31_2018-10-09/escherichia_coli_pubmlst_k31_m023_2018-10-09.jld' -1 input_t1.fq")
    os.system("mv mentalist_out.byvote " + args.mlstsevenloci)
    sequence_typing = openFileAsTable(args.mlstsevenloci)
    sequence_qc = openFileAsTable("mentalist_out.coverage.txt")
    log.write("\n\nSequence Typer\n==============\n")
    log.write(os.popen("mentalist -v | grep MentaLiST").read())
    log.write("\n")
    log.write(os.popen("cat " +  TOOL_DIR + "/data/SequenceTyping_db.txt").read())
    if args.shigatoxintyping:
        # SHIGATOXIN TYPER
        os.system("echo 'sseqid\tpident\tlength\tpositive' > mmseqs_shigatoxin_fct")
        os.system("mmseqs easy-search --search-type 3 --format-output query,target,pident,alnlen,tlen contigs.fa " + TOOL_DIR + "/data/stxDB mmseqs_STX /tmp")
        os.system(TOOL_DIR + "/scripts/extract_shigatoxin.sh")
        os.system("cat mmseqs_shigatoxin_fc >> mmseqs_shigatoxin_fct")
        shigatoxin_typing = openFileAsTable("mmseqs_shigatoxin_fc")
        log.write("\n\nShigatoxin Typer\n==============\n")
        log.write(os.popen("cat " +  TOOL_DIR + "/data/ShigatoxinTyping_db.txt").read())
    if args.serotyping:        
        # SEROTYPER
        os.system("echo 'sseqid\tpident\tlength\tpositive' > mmseqs_OH_fc")
        # SEROTYPER O
        os.system("mmseqs easy-search --search-type 3 --format-output target,pident,alnlen,tlen contigs.fa " + TOOL_DIR + "/data/O_typeDB mmseqs_O /tmp")
        os.system("awk -F '\t' '($3>800 && $4>800) { print $1 FS $2 FS $3 FS $4 }' mmseqs_O | sort -nrk 2 -nrk 3 > mmseqs_O_fc")
        sero_typing_o = openFileAsTable("mmseqs_O_fc")
        os.system("cat mmseqs_O_fc >> mmseqs_OH_fc")
        # SEROTYPER H
        os.system("mmseqs easy-search --search-type 3 --format-output target,pident,alnlen,tlen contigs.fa " + TOOL_DIR + "/data/H_typeDB mmseqs_H /tmp")
        os.system("awk -F '\t' '($3>800 && $4>800) { print $1 FS $2 FS $3 FS $4 }' mmseqs_H | sort -nrk 2 -nrk 3 > mmseqs_H_fc")
        sero_typing_h = openFileAsTable("mmseqs_H_fc")
        os.system("cat mmseqs_H_fc >> mmseqs_OH_fc")    
        if os.stat('mmseqs_O_fc').st_size == 0 and os.stat('mmseqs_H_fc').st_size == 0:
            os.system("echo '-\t-\t-\t-' >> mmseqs_OH_fc")
        log.write("\n\nSero Typer\n==============\n")
        log.write(os.popen("cat " +  TOOL_DIR + "/data/SeroTyping_db.txt").read())
    if args.amrtyping:
        # AMRGENES
        os.system("amrfinder --threads 4 --database " + TOOL_DIR + "/data/amrfinder -n contigs.fa -O Escherichia -o " + args.amr)
        log.write("\n\nAMR Typer\n==============\nAMRFinderPlus ")
        log.write(os.popen("amrfinder --version").read())
        log.write("\ndatabase version: ")
        log.write(os.popen("cat " + TOOL_DIR + "/data/amrfinder/version.txt").read())
    # REPORT
    try:
        report = open(args.output, 'w')
        # write head html
        insertFile(TOOL_DIR + "/report_head.html", report)
        report.write("<td><h1>EURL VTEC WGS PT</h1><h2>Report for %s</h2>%s</td>" % (args.input1_name, datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC"))) 
        insertFile(TOOL_DIR + "/report_head2.html", report)
        # write results
        report.write("<h3>Summary</h3>\n")
        if args.serotyping:
            report.write("<p>Serotype: ")
            if len(sero_typing_o) == 0:
                report.write("O?")
            else:
                report.write("%s" % sero_typing_o[0][0][sero_typing_o[0][0].rfind("O"):])
            if len(sero_typing_h) == 0:
                report.write(":H?")
            else:
                report.write(":%s" % sero_typing_h[0][0][sero_typing_h[0][0].rfind("H"):])
            report.write("</p>\n")
        report.write("<p>Sequence type: ")
        if len(sequence_typing) < 2:
            report.write("Sequence typing failed")
        elif sequence_typing[1][1] == "failed":
            report.write("Sequence typing failed")
        else:
            report.write("ST%s" % sequence_typing[1][8])
        report.write("</p>\n")
        if args.virulotyping:
            os.system("sort " + args.virulotyper  + " | awk '/eae_|stx1._|stx2._|ehxa_/ && $2>50 && !seen[substr($1, 1, index($1, \"_\")-2)]++ { printf(\"%s%s\",sep,substr($1, 1, index($1, \"_\")-1));sep=\", \" }END{print \"\"}' > virulotyper_rep")
            for line in fileinput.input("virulotyper_rep", inplace=True):
                print(line.replace("1a", "1"),)
            for line in fileinput.input("virulotyper_rep", inplace=True):
                print(line.replace("2a", "2"),)
            for line in fileinput.input("virulotyper_rep", inplace=True):
                print(line.replace("1b", "1"),)
            for line in fileinput.input("virulotyper_rep", inplace=True):
                print(line.replace("2b", "2"),)
            report.write("<p>Virulotypes: ")
            insertFile("virulotyper_rep", report)
            report.write("</p>\n")
        if args.shigatoxintyping:
            report.write("<p>Stx Subtypes: ")
            if len(shigatoxin_typing) == 0:
                report.write("No subtype match found")
            else:
                shigatoxin_subtype = ""
                shigatoxin_types = openFileAsTable(TOOL_DIR + "/data/stx_subtypes")
                for subtype in shigatoxin_typing:
                    blast_pident_100 = float(subtype[1]) == 100
                    if (blast_pident_100):
                        for item in shigatoxin_types:
                            if item[0] == subtype[0]:
                                shigatoxin_subtype = item[1] + " " + shigatoxin_subtype
                if len(shigatoxin_subtype) == 0:
                    shigatoxin_subtype = "No complete subtype match found"
                report.write("%s" % shigatoxin_subtype)
            report.write("</p>\n")
        # Quality Check
        disclaimer = False
        sequence_qc_cov = 0
        for x in range(7):
            if float(sequence_qc[x+1][2]) < 1:
                disclaimer = True
            sequence_qc_cov = sequence_qc_cov + float(sequence_qc[x+1][3])
        if sequence_qc_cov < 210:
            disclaimer = True
        if disclaimer:
            report.write("<p style='font-weight:bold;color:red'>Disclaimer: The data analysed do not fulfill minimum quality parameters, please consider repeating the sequencing.</p>\n")
        report.write("<hr/><h3>Raw data quality check</h3>\n")
        if args.input2:
            report.write("<p>FASTQC result forward: <a href='%s/datasets/%s/display/?preview=True'>Webpage</a></p>\n" % (BASE_URL, args.html1_id))
            report.write("<p>FASTQC result reverse: <a href='%s/datasets/%s/display/?preview=True'>Webpage</a></p>\n" % (BASE_URL, args.html2_id))
        else:
            report.write("<p>FASTQC result: <a href='%s/datasets/%s/display/?preview=True'>Webpage</a></p>\n" % (BASE_URL, args.html1_id))
        if args.serotyping:
            report.write("<br/><hr/><h3>Serotyping</h3>\n")
            insertFileAsTable("mmseqs_OH_fc", report, True)
        report.write("<br/><hr/><h3>Multi Locus Sequence Typing</h3>\n")
        if len(sequence_typing) > 1:
            insertTable(sequence_typing, report, True)
        if args.virulotyping:
            report.write("<br/><hr/><h3>Virulotyping</h3>\n")
            report.write("<p>This table is filtered for results with >90%% gene coverage, unfiltered results can be found <a href='%s/datasets/%s/display/?preview=True'>here</a></p>\n" % (BASE_URL, args.virulotyper_id))
            insertFileAsTable("pathotyper_rep_tab", report, True, "table table-cross")
        if args.shigatoxintyping:
            report.write("<br/><hr/><h3>Shiga toxin typing</h3>\n")
            insertFileAsTable("mmseqs_shigatoxin_fct", report, True)
        if args.amrtyping:
            report.write("<br/><hr/><h3>AMR typing</h3>\n")
            report.write("<p>AMR result: <a href='%s/datasets/%s/display/?preview=True'>Webpage</a></p>\n" % (BASE_URL, args.amr_id))
        # write tail html
        insertFile(TOOL_DIR + "/report_tail.html", report)
    finally:
        report.close()

if __name__ == "__main__":
    __main__()
