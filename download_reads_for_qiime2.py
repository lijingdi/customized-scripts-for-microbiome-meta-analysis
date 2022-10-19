import re
import sys
import os
import subprocess
import logging
import argparse
from baseops import file


project = sys.argv[1]
study = sys.argv[2]
path = "/data/zool-temp-micro/pemb6000/16S/" + study
meta = open(path + "/" + project + ".tsv", "r")

next(meta)
for line in meta:
	line = line.strip()
	col = line.split("\t")
	links = col[20]
	run = col[5]
	lib_name_raw = col[-2]
	lib_name = "-".join(lib_name_raw.split(".")) + "_L001"
	fq_name = run + "_" + lib_name
	if ";" in links:
		link1 = links.split(";")[0]
		link2 = links.split(";")[1]
		fq1 = path + "/raw/" + fq_name + "_R1_001.fastq.gz"
		fq2 = path + "/raw/" + fq_name + "_R2_001.fastq.gz"
		try:
			with open(path + "/" + project + "_download.log", "a") as fout:
				if not file.exists(fq1):
					subprocess.run("curl -o " + fq1 + " " + link1, check=True, shell=True, stdout=fout, stderr=fout)
				if not file.exists(fq2):
					subprocess.run("curl -o " + fq2 + " " + link2, check=True, shell=True, stdout=fout, stderr=fout)			
		except subprocess.CalledProcessError:
			logging.info("download fails, check logfile {}". format(fout))
	else:
		try:
			with open(path + "/" + project + "_download.log", "a") as fout:
				fq = path + "/raw/" + fq_name + "_R1_001.fastq.gz"
				if not file.exists(fq):
					subprocess.run("curl -o " + fq + " " + links, check=True, shell=True, stdout=fout, stderr=fout)
		except subprocess.CalledProcessError:
			logging.info("download fails, check logfile {}". format(fout))


