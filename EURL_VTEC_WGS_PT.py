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
import shutil
import subprocess
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
    parser.add_argument('--contigs', dest='contigs', help='Assembly contigs')
    parser.add_argument('--quast', dest='quast', help='Quast report')
    parser.add_argument('--log', dest='logfile', help='log file')
    parser.add_argument('--virulotyper', dest='virulotyper', help='Virulotyping Mapping reads')
    parser.add_argument('--virulotyper_id', dest='virulotyper_id', help='Virulotyping Mapping reads id')
    parser.add_argument('--stx', dest='stx', help='Shiga toxin')
    parser.add_argument('--mlstsevenloci', dest='mlstsevenloci', help='Multi Locus Alleles table')
    parser.add_argument('--amr', dest='amrgenes', help='AMR genes')
    parser.add_argument('--amr_id', dest='amr_id', help='AMR file id')
    parser.add_argument('--antigen_O', dest='antigen_O', help='Antigen for O')
    parser.add_argument('--antigen_H', dest='antigen_H', help='Antigen for H')
    parser.add_argument('--output', dest='output', help='output report html file')
    args = parser.parse_args()

    log = open(args.logfile, 'w')
    log.write("EURL VTEC WGS PT v3.2\n\nTool versions\n=============\n")
    os.system("ln -s " + os.popen("$(dirname $(readlink -e $(which trimmomatic)))").read().strip() + "/trimmomatic.jar trimmomatic.jar")
    # FASTQC
    subprocess.call("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input1 + " -d " + args.html1_path + " -o " + args.html1 + " -t " + args.text1 + " -f " + args.input1_ext + " -j " + args.input1_name + " -e " + "fastqc", shell=True)
    log.write(os.popen("fastqc -v").read())
    if args.input2:
        # FASTQC
        subprocess.call("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input2 + " -d " + args.html2_path + " -o " + args.html2 + " -t " + args.text2 + " -f " + args.input2_ext + " -j " + args.input2_name + " -e " + "fastqc", shell=True)
        # TRIMMING
        subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar PE -threads ${GALAXY_SLOTS:-6} -phred33 " + args.input1 + " " + args.input2 + " trimmed1.fq trimmed1unpaired trimmed2.fq trimmed2unpaired SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36", shell=True)
        log.write("\nTrimmomatic v0.39\n")
        log.write("parameters: phred33 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36\n\n")
        # ASSEMBLY
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl spades_contigs spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} --pe1-ff --pe1-1 trimmed1.fq --pe1-2 trimmed2.fq", shell=True)
        subprocess.call("perl " + TOOL_DIR + "/scripts/filter_spades_repeats.pl -i spades_contigs -t spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o output_with_repeats -u output_without_repeats -n repeat_sequences_only -e 5000 -f discarded_sequences -s summary", shell=True)
        shutil.move("output_without_repeats", args.contigs)
        log.write(os.popen("spades.py -v").read())
        log.write("parameters: --isolate, pe1-ff, pe1-1, pe1-2 filter_repeats\n\n")
    else:
        # TRIMMING
        subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 " + args.input1 + " trimmed1.fq SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55", shell=True)
        log.write("\nTrimmomatic v0.39\n")
        log.write("parameters: phred33 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55\n\n")
        # ASSEMBLY
        subprocess.call("perl " + TOOL_DIR + "/scripts/spades.pl spades_contigs spades_contig_stats spades_scaffolds spades_scaffold_stats spades_log NODE spades.py --disable-gzip-output --isolate -t ${GALAXY_SLOTS:-16} --iontorrent -s trimmed1.fq", shell=True)
        subprocess.call("perl " + TOOL_DIR + "/scripts/filter_spades_repeats.pl -i spades_contigs -t spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o output_with_repeats -u output_without_repeats -n repeat_sequences_only -e 5000 -f discarded_sequences -s summary", shell=True)
        shutil.move("output_without_repeats", args.contigs)
        log.write(os.popen("spades.py -v").read())
        log.write("parameters: --isolate, --iontorrent filter_repeats\n\n")
    # QUAST
    subprocess.call("quast --threads 4 -o outputdir --est-ref-size 5000000 --min-contig 500 -l  '" + args.input1_name + "' --contig-thresholds 0,1000 " + args.contigs, shell=True)
    shutil.move("outputdir/report.tsv", args.quast)
    if args.virulotyping:
        # VIRULOTYPER
        if args.input2:
            subprocess.call("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f " + args.input1 + " " + args.input2 + " -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        else:
            subprocess.call("perl " + TOOL_DIR + "/scripts/patho_typing.pl 'python " + TOOL_DIR + "/scripts/patho_typing.py -s Escherichia coli -f " + args.input1 + " -o output_dir -j 4 --minGeneCoverage 90 --minGeneIdentity 90 --minGeneDepth 15'", shell=True)
        subprocess.call("(head -n 1 pathotyper_rep_tot_tab && tail -n +2 pathotyper_rep_tot_tab | sort -k 2rn) > " + args.virulotyper, shell=True)
        log.write("\n\nViruloTyper\n===========\npatho_typing v1.0\n")
        log.write("parameters: minGeneCoverage=90, minGeneIdentity=90, minGeneDepth=15\n\n")
        log.write(os.popen("cat " +  TOOL_DIR + "/data/ViruloTyping_db.txt").read())
    # SEQUENCETYPER
    subprocess.call("mlst --legacy --scheme ecoli " + args.contigs + " | cut -f3,4,5,6,7,8,9,10 > " + args.mlstsevenloci, shell=True)
    sequence_typing = openFileAsTable(args.mlstsevenloci)
    log.write("\n\nSequence Typer\n==============\n")
    log.write(os.popen("mlst -v").read())
    log.write("\n")
    log.write(os.popen("cat " +  TOOL_DIR + "/data/SequenceTyping_db.txt").read())
    if args.shigatoxintyping:
        # SHIGATOXIN TYPER
        if args.input2:
            # CONSENSUS
            subprocess.call("sh " + TOOL_DIR + "/scripts/stx_subtype_pe.sh " + TOOL_DIR + " trimmed1.fq trimmed2.fq " + args.contigs, shell=True)
        else:
            # CONSENSUS
            subprocess.call("sh " + TOOL_DIR + "/scripts/stx_subtype_se.sh " + TOOL_DIR + " trimmed1.fq " + args.contigs, shell=True)
        # SHIGATOXIN SEQUENCE SEARCH
        subprocess.call("sh " + TOOL_DIR + "/scripts/stx_subtype_fa.sh " + TOOL_DIR + " stx.fasta", shell=True)
        subprocess.call("echo 'sseqid\tpident\tlength\tpositive' > shigatoxin_fct", shell=True)
        subprocess.call("cat shigatoxin_fc >> shigatoxin_fct", shell=True)
        shutil.move("shigatoxin_fct", args.stx)
        shigatoxin_typing = openFileAsTable("shigatoxin_fc")
        log.write("\n\nShigatoxin Typer v2.0\n==============\n")
        log.write(os.popen("cat " +  TOOL_DIR + "/data/ShigatoxinTyping_db.txt").read())
    if args.serotyping:        
        # SEROTYPER
        subprocess.call("echo 'sseqid\tpident\tlength\tpositive' > serogroup_OH_fcd", shell=True)
        if args.input2:
            subprocess.call("sh " + TOOL_DIR + "/scripts/serotype.sh " + TOOL_DIR + " y " + args.input1 + " " + args.input2 + " " + args.contigs, shell=True)
        else:
            subprocess.call("sh " + TOOL_DIR + "/scripts/serotype.sh " + TOOL_DIR + " n " + args.input1 + " xxx " + args.contigs, shell=True)
        # SEROTYPER O
        subprocess.call("awk -F '\t' '$4>800 { print $2 FS $3 FS $4 FS $16 }' serogroup_O | sort -nrk 2 -nrk 3 > serogroup_O_fc", shell=True)
        subprocess.call("awk -F , '!seen[$0]++' serogroup_O_fc > serogroup_O_fcd", shell=True)
        sero_typing_o = openFileAsTable("serogroup_O_fcd")
        subprocess.call("cat serogroup_O_fcd >> serogroup_OH_fcd", shell=True)
        shutil.move("serogroup_O_fcd", args.antigen_O)
        # SEROTYPER H
        subprocess.call("awk -F '\t' '$4>800 { print $2 FS $3 FS $4 FS $16 }' serogroup_H | sort -nrk 2 -nrk 3 > serogroup_H_fc", shell=True)
        subprocess.call("awk -F , '!seen[$0]++' serogroup_H_fc > serogroup_H_fcd", shell=True)
        sero_typing_h = openFileAsTable("serogroup_H_fcd")
        subprocess.call("cat serogroup_H_fcd >> serogroup_OH_fcd", shell=True)
        shutil.move("serogroup_H_fcd", args.antigen_H)
        if os.stat(args.antigen_O).st_size == 0 and os.stat(args.antigen_H).st_size == 0:
            subprocess.call("echo '-\t-\t-\t-' >> serogroup_OH_fcd", shell=True)
        log.write("\n\nSero Typer\n==============\n")
        log.write(os.popen("cat " +  TOOL_DIR + "/data/SeroTyping_db.txt").read())
    if args.amrtyping:
        # AMRGENES
        subprocess.call("amrfinder --threads 4 --database " + TOOL_DIR + "/data/amrfinder -n " + args.contigs + " -O Escherichia -o " + args.amrgenes, shell=True)
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
        elif sequence_typing[1][1] == "-":
            report.write("Sequence typing failed")
        else:
            report.write("ST%s" % sequence_typing[1][0])
        report.write("</p>\n")
        if args.virulotyping:
            subprocess.call("sort " + args.virulotyper  + " | awk '/eae_|stx1._|stx2._|ehxa_/ && $2>50 && !seen[substr($1, 1, index($1, \"_\")-2)]++ { printf(\"%s%s\",sep,substr($1, 1, index($1, \"_\")-1));sep=\", \" }END{print \"\"}' > virulotyper_rep", shell=True)
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
                str_shigatoxin_subtype = "No subtype match found"
            else:
                # get corresponding subtypes
                str_shigatoxin_subtype = ""
                shigatoxin_subtypes = []
                shigatoxin_subtypes_raw = []
                shigatoxin_types = openFileAsTable(TOOL_DIR + "/data/stx_subtypes")
                for subtype in shigatoxin_typing:
                    blast_pident_100 = float(subtype[1]) == 100
                    if (blast_pident_100):
                        for item in shigatoxin_types:
                            if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                                shigatoxin_subtypes.append(item[1])
                                shigatoxin_subtypes_raw.append(item[1])
               # partial matches
                for subtype in shigatoxin_typing:
                    for item in shigatoxin_types:
                        if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                            if item[1][0:4] == "stx1":
                                shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                                shigatoxin_subtypes_raw.append(item[1])
                            if item[1][0:4] == "stx2":
                                shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                                shigatoxin_subtypes_raw.append(item[1])
                shigatoxin_subtypes.sort()
                str_shigatoxin_subtype = " ".join(shigatoxin_subtypes)
            report.write("%s" % str_shigatoxin_subtype)
            report.write("</p>\n")
        # Quality Check
        disclaimer = False
        if any("-" in s for s in sequence_typing) or any("?" in s for s in sequence_typing):
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
            insertFileAsTable("serogroup_OH_fcd", report, True)
        report.write("<br/><hr/><h3>Multi Locus Sequence Typing</h3>\n")
        if len(sequence_typing) > 1:
            insertTable(sequence_typing, report, True)
        if args.virulotyping:
            report.write("<br/><hr/><h3>Virulotyping</h3>\n")
            report.write("<p>This table is filtered for results with >90%% gene coverage, unfiltered results can be found <a href='%s/datasets/%s/display/?preview=True'>here</a></p>\n" % (BASE_URL, args.virulotyper_id))
            insertFileAsTable("pathotyper_rep_tab", report, True, "table table-cross")
        if args.shigatoxintyping:
            report.write("<br/><hr/><h3>Shiga toxin typing</h3>\n")
            insertFileAsTable(args.stx, report, True)
        if args.amrtyping:
            report.write("<br/><hr/><h3>AMR typing</h3>\n")
            report.write("<p>AMR result: <a href='%s/datasets/%s/display/?preview=True'>Webpage</a></p>\n" % (BASE_URL, args.amr_id))
        # write tail html
        insertFile(TOOL_DIR + "/report_tail.html", report)
    finally:
        report.close()

if __name__ == "__main__":
    __main__()
