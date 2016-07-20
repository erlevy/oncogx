#!/usr/bin/python
import os
import subprocess
import csv
import xml.etree.ElementTree as ET

output = "/Users/Eric/tcga/results/"
XML = "/Users/Eric/tcga/XML/"
measurements = "/Users/Eric/tcga/clinical/"

# list all files in this directory
files = os.listdir(measurements)

data_table = list()
annai_table = list()
data_table.append(["Patient UUID", "Sample Barcode", "Center Name", "Slide ID", "Aliquot ID", "Cancer Type", "Percent Lymphocytes", "Reference Name", "Exome ID", "File Size", "Checksum", "Upload Date"])
# list of lists with each component having"
# patient UUID, sample barcode, lymphocyte percent,(from biospecimen slide)
# tumor type
# exome UUID
#nationwidechildrens.org_biospecimen_slide_blca
# RNA-seq UUID
for name in files:
	# get patient ID, sample barcode, cancer type, and percent lymphocytes from biospecimen data
	type = name[42:]
	type = type[:-4]
	id == ""
	if (name[0] == "." or name == "nationwidechildrens.org_biospecimen_slide_fppp.txt" or name == "nationwidechildrens.org_biospecimen_slide_meso.txt" or name == "nationwidechildrens.org_biospecimen_slide_uvm.txt" or name == "nationwidechildrens.org_biospecimen_slide_sarc.txt" or name == "nationwidechildrens.org_biospecimen_slide_chol.txt" or name == "nationwidechildrens.org_biospecimen_slide_tgct.txt" or name == "nationwidechildrens.org_biospecimen_slide_thym.txt"):
		continue
	with open(measurements + name, "rb") as f:
		reader = csv.reader(f, delimiter="\t")
		for row in reader:
			if (row[0] != "bcr_patient_uuid" and row[0] != "" and row[6] != '[Not Available]'):
				id = row[0]
				sample = row[1]
				slide = row[2]
				if (sample[-3:-1] != "01"):
					continue
				# get XLM files for each patient
				participant_id = "participant_id=" + id
#				cgquery_output = subprocess.check_output(['cgquery', participant_id, '-o', XML + id + '.xml'])
				input_xml = XML + id + '.xml'
				tree = ET.parse(input_xml)
				root = tree.getroot()
				wxs_id = ""
				sample_id = ""
				reference_name = ""
				data_type = ""
				file_size = ""
				checksum = ""
				center_name = ""
				aliquot_id = ""
				cancer_type = ""
				upload_date = ""
				for result in root.iter ("Result"):
					analysis_id = ""
					is_wxs = 0
					is_mapped = 0
					is_tumor = 0
					is_tcga = 0
					is_matched = 0
					is_illumina = 0
					for child in result:
						if (child.tag=="analysis_id"):
							analysis_id = child.text
						if (child.tag=="library_strategy"):
							if (child.text=="WXS"):
								is_wxs = 1
								data_type = child.text
						if (child.tag=="center_name"):
							center_name = child.text
						if (child.tag=="aliquot_id"):
							aliquot_id = child.text
						if (child.tag=="disease_abbr"):
							cancer_type = child.text
						if (child.tag=="upload_date"):
							upload_date = child.text
						if (child.tag=="refassem_short_name"):
							if (child.text!="unaligned"):
								is_mapped = 1
								reference_name = child.text
						if (child.tag=="sample_type"):
							if (child.text=="01"):
								is_tumor = 1
						if (child.tag=="study"):
							if (child.text=="phs000178"):
								is_tcga = 1
						if (child.tag=="platform"):
							if (child.text=="ILLUMINA"):
								is_illumina = 1
						if (child.tag=="files"):
							for file in child:
								is_bam = 0
								for file_info in file:
									if (file_info.tag=="filename"):
										if (file_info.text[-4:]==".bam"):
											is_bam = 1
									if (file_info.tag=="filesize" and is_bam==1):
										file_size = file_info.text
									if (file_info.tag=="checksum" and is_bam==1):
										checksum = file_info.text
						if (child.tag=="legacy_sample_id"):
							if (sample in child.text):
								is_matched = 1
								sample_id = child.text
						if (child.tag=="legacy_sample_id"):
							sample_id = child.text
						if (is_wxs==1 and is_mapped==1 and is_tumor==1 and is_tcga==1 and is_matched==1 and is_illumina==1):
							wxs_id = analysis_id
							ids_row = [row[0], sample_id, center_name, slide, aliquot_id, cancer_type, row[6], reference_name, wxs_id, round(float(file_size)/1000000000, 2), checksum, upload_date]
							if (ids_row in data_table or slide[-3:-2] == "M"):
								pass
							else:
								data_table.append(ids_row)
							is_wxs=0
							is_mapped=0
							is_tumor=0
							is_tcga=0
							is_matched=0
							is_illumina=0
							analysis_id = ""

	f.close()
			
with open(output + "lymphocyte_data_table.csv", "wb") as g:
	writer = csv.writer(g)
	writer.writerows(data_table)

