from genomeload import load_genome
from gtfline import gtf_line
import pandas as pd
import numpy as np
import gzip
import sys

print("Start loading Genome")
#file path for genome sequence
genome = load_genome("../pri_hg38.fa")

print("Finish loading Genome")

print("Start loading annotated gtf file")
rna_set = {}

# File path for the annotated GTF file
annotation_file = "../gencode.v38.annotation.gtf"

# Initialize rna_set with keys for both strands of each chromosome in genome
for key in genome:
    rna_set[key + "+"] = {}
    rna_set[key + "-"] = {}

# Open the annotation file, supporting gzipped files
try:
    if annotation_file.endswith(".gz"):
        af = gzip.open(annotation_file, 'rt')  # 'rt' mode for text
    else:
        af = open(annotation_file, 'r')
except IOError as e:
    print(f"The annotation file could not be opened: {e}")
    sys.exit(1)

tc = 0

with af:
    for content in af:
        fields = content.strip().split("\t")
        if len(fields) < 9:  # Ensure there are at least 9 fields
            continue

        gtf_list = gtf_line(content.strip())
        chrom_strand_key = gtf_list.chrom + gtf_list.strand

        if gtf_list.type == "gene":
            gene_name = gtf_list.attributes["gene_name"]
            rna_set[chrom_strand_key][gene_name] = {}

        elif gtf_list.type == "transcript":
            gene_name = gtf_list.attributes["gene_name"]
            transcript_id = gtf_list.attributes["transcript_id"]
            rna_set[chrom_strand_key][gene_name][transcript_id] = {
                "TSS": gtf_list.start,
                "TES": gtf_list.end,
                "exons": [],
                "start_codon":[],
                "CDS": []
            }
            tc += 1
            if tc % 10000 == 0:
                print(f"Loaded {tc} transcripts")

        elif gtf_list.type == "exon":
            gene_name = gtf_list.attributes["gene_name"]
            transcript_id = gtf_list.attributes["transcript_id"]
            rna_set[chrom_strand_key][gene_name][transcript_id]["exons"].append([gtf_list.start, gtf_list.end])

        elif gtf_list.type == "start_codon":
            gene_name = gtf_list.attributes["gene_name"]
            transcript_id = gtf_list.attributes["transcript_id"]
            if not rna_set[chrom_strand_key][gene_name][transcript_id]["start_codon"]:
                rna_set[chrom_strand_key][gene_name][transcript_id]["start_codon"] = [gtf_list.start, gtf_list.end]
            else:
                current_start, current_end = rna_set[chrom_strand_key][gene_name][transcript_id]["start_codon"]
                updated_start = min(current_start, gtf_list.start)
                updated_end = max(current_end, gtf_list.end)
                rna_set[chrom_strand_key][gene_name][transcript_id]["start_codon"] = [updated_start, updated_end]

        elif gtf_list.type == "CDS":
            gene_name = gtf_list.attributes["gene_name"]
            transcript_id = gtf_list.attributes["transcript_id"]
            cds_list = rna_set[chrom_strand_key][gene_name][transcript_id]["CDS"]
            if cds_list:
                cds_list[0] = min(cds_list[0], gtf_list.start)
                cds_list[1] = max(cds_list[1], gtf_list.end)
            else:
                rna_set[chrom_strand_key][gene_name][transcript_id]["CDS"] = [gtf_list.start, gtf_list.end]
print("Load finished")

print("Start loading classify file")

#file path for classify file
cf = "./293T_R1_classification.txt"
df = pd.read_csv(cf, delimiter="\t")

# Set 'isoform' as the index of the DataFrame
df.set_index('isoform', inplace=True)
df['fl_assoc'] = pd.to_numeric(df['fl_assoc'], errors='coerce').astype('Int64')


print("Start loading new assembled GFF file")
rna_set2 = {}

#file path for new assembled gtf file
annotation_file = "../293T_R1.gff"

# Initialize rna_set with keys for both strands of each chromosome in genome
for key in genome:
    rna_set2[key + "+"] = {}
    rna_set2[key + "-"] = {}

# Open the annotation file, supporting gzipped files
try:
    if annotation_file.endswith(".gz"):
        af = gzip.open(annotation_file, 'rt')  # 'rt' mode for text
    else:
        af = open(annotation_file, 'r')
except IOError as e:
    print(f"The annotation file could not be opened: {e}")
    sys.exit(1)

content = af.readline()
tc = 0

with af:
    for content in af:
        fields = content.strip().split("\t")
        if len(fields) < 9:  # Ensure there are at least 9 fields
            continue

        gtf_list = gtf_line(content.strip())
        chrom_strand_key = gtf_list.chrom + gtf_list.strand

        if gtf_list.type == "transcript":
            
            gene_id = gtf_list.attributes["gene_id"]
            transcript_id = gtf_list.attributes["transcript_id"]
            if gene_id not in rna_set2[chrom_strand_key]:
                rna_set2[chrom_strand_key][gene_id] = {}
            
            rna_set2[chrom_strand_key][gene_id][transcript_id] = {
                "TSS": gtf_list.start,
                "TES": gtf_list.end,
                "exons": []
            }
            tc += 1
            if tc % 10000 == 0:
                print(f"Loaded {tc} transcripts")
                
        elif gtf_list.type == "exon":
            rna_set2[chrom_strand_key][gene_id][transcript_id]["exons"].append([gtf_list.start, gtf_list.end])

print("Load new assembled GFF file finished")

# Specify the output file path
output_file_path = "./293.xlsx"
"""
with pd.ExcelWriter(output_file_path, mode='w', engine='openpyxl') as writer:
    empty_df = pd.DataFrame()
    empty_df.to_excel(writer, sheet_name='Placeholder')
"""
print("Start loading target anisoform information")

#file path for target gene file
tf = "./target_R1.csv"
df2 = pd.read_csv(tf, delimiter=',')

# Set 'Gene Name' as the index of the DataFrame
df2.set_index('Gene Name', inplace=True)

# Initialize empty DataFrames for each output CSV
df_fl_assoc_counts_by_level = pd.DataFrame()
df_fl_assoc_percentages_by_level = pd.DataFrame()
df_isoform_counts_by_level = pd.DataFrame()
df_isoform_percentages_by_level = pd.DataFrame()

# Initialize empty DataFrames for coding statistics
df_fl_assoc_counts_by_coding = pd.DataFrame()
df_fl_assoc_percentages_by_coding = pd.DataFrame()
df_isoform_counts_by_coding = pd.DataFrame()
df_isoform_percentages_by_coding = pd.DataFrame()

# Define all possible coding types
coding_types = ["CDS", "anisoform", "non_coding", "no Stop"]

for gene_name, row in df2.iterrows():
    if gene_name not in df["associated_gene"].values:
        continue

    # Extract ateORF coordinates information
    aORF = row["anisoform coordinates(hg38)"]
    coordinates = aORF.split(";")
    chkey = coordinates[0].split(":")[0] + row["Strand"]
    start_set, end_set = set(), set()
    for coord in coordinates:
        chrom, positions = coord.split(":")
        start_pos, end_pos = positions.split("-")
        # Convert start and end positions to integers and add to sets
        start_set.add(int(start_pos.replace(',', '')))
        end_set.add(int(end_pos.replace(',', '')))
    
    sset0, sset1 = set(),set()
    for t,inf in rna_set[chkey][gene_name].items():
        if inf["start_codon"]:
            sset0.add(inf["start_codon"][0])
            sset1.add(inf["start_codon"][1])

    if "-" in chkey:
        cdssp = sset1
        asp = end_set
    else:
        cdssp = sset0
        asp = start_set

    # Create tuples with the source information
    cdssp_tuples = [(x, 'CDS') for x in cdssp]
    asp_tuples = [(x, 'anisoform') for x in asp]

    # Merge the two lists of tuples
    #merged_splist = cdssp_tuples
    merged_splist = cdssp_tuples + asp_tuples

    # Sort the merged list by the numeric value
    sorted_splist = sorted(merged_splist, key=lambda x: x[0], reverse=("-" in chkey))
    
    aSeq = row["Sequence"]
    lines = df[df["associated_gene"] == gene_name]
    lines = lines.sort_values(by='fl_assoc', ascending=False)
    lines['sort_order'] = range(1, len(lines) + 1)
    total_fl_assoc = lines['fl_assoc'].sum()
    lines['percentage'] = lines['fl_assoc'] / total_fl_assoc * 100

    conditions = [
    (lines['sort_order'] <= 10) & (lines['percentage'] >= 1),
    (lines['sort_order'] > 10) & (lines['percentage'] >= 1),
    (lines['sort_order'] <= 10) & (lines['percentage'] < 1)
    ]
    
    choices = [
    'Significant isoform',
    'Isoform above 1%',
    'Top 10 isoform'
    ]
    
    lines['Evidence level'] = np.select(conditions, choices, default='Detectable isoform')

    # Calculate the sum of fl_assoc for each Evidence level
    fl_assoc_counts = lines.groupby('Evidence level')['fl_assoc'].sum()

    # Calculate the percentage of the total fl_assoc
    fl_assoc_percentages = fl_assoc_counts / total_fl_assoc * 100

    # Count the number of isoforms in each Evidence level
    isoform_counts = lines['Evidence level'].value_counts()

    # Calculate the percentage of isoform counts
    total_isoforms = len(lines)
    isoform_percentages = isoform_counts / total_isoforms * 100

    # Convert results to DataFrame and add gene_name as a column for context
    fl_assoc_counts_df = pd.DataFrame(fl_assoc_counts).reset_index()
    fl_assoc_counts_df['gene_name'] = gene_name

    fl_assoc_percentages_df = pd.DataFrame(fl_assoc_percentages).reset_index()
    fl_assoc_percentages_df['gene_name'] = gene_name

    isoform_counts_df = pd.DataFrame(isoform_counts).reset_index()
    isoform_counts_df['gene_name'] = gene_name

    isoform_percentages_df = pd.DataFrame(isoform_percentages).reset_index()
    isoform_percentages_df['gene_name'] = gene_name

    # Append each DataFrame to the corresponding summary DataFrame
    df_fl_assoc_counts_by_level = pd.concat([df_fl_assoc_counts_by_level, fl_assoc_counts_df])
    df_fl_assoc_percentages_by_level = pd.concat([df_fl_assoc_percentages_by_level, fl_assoc_percentages_df])
    df_isoform_counts_by_level = pd.concat([df_isoform_counts_by_level, isoform_counts_df])
    df_isoform_percentages_by_level = pd.concat([df_isoform_percentages_by_level, isoform_percentages_df])
    
    Gene_id = ".".join(lines.index[0].split(".")[:-1])
    anno_trans = rna_set[chkey][gene_name]
    assem_trans = rna_set2[chkey][Gene_id]
    
    for iso, line in lines.iterrows():
        asso_tran = line.loc["associated_transcript"]
        exons =  assem_trans[iso]["exons"]
        exons.sort()
        POSLIST, seqLIST = [], ""
        for exon in exons:
            POSLIST += range(exon[0], exon[1] + 1)
            seqLIST += genome[chkey[:-1]][exon[0]-1:exon[1]]
        if chkey[-1] == "-":
            POSLIST = POSLIST[::-1]
            seqLIST = seqLIST.reverse_complement()
        for p,source in sorted_splist:
            if p not in POSLIST:
                pass
            else:
                POSslice = POSLIST[POSLIST.index(p):]
                SEQslice = seqLIST[POSLIST.index(p):]
                SEQslice = SEQslice[:(len(SEQslice) // 3) * 3] 
                pep = SEQslice.translate()
                if "*" in pep:
                    lines.at[iso, "coding"] = source
                    break
                else:
                    lines.at[iso, "coding"] = "no Stop"
    with pd.ExcelWriter(output_file_path, mode='a', engine='openpyxl') as writer:
        lines.to_excel(writer, sheet_name=gene_name, index=True, na_rep='NaN')
    print(f"Results for {gene_name} have been saved.")

    # Calculate statistics for coding
    coding_counts = lines.groupby('coding')['fl_assoc'].sum()
    coding_percentages = coding_counts / total_fl_assoc * 100

    isoform_counts = lines['coding'].value_counts()
    isoform_percentages = isoform_counts / total_isoforms * 100

    coding_counts = coding_counts.reindex(coding_types, fill_value=0)
    coding_percentages = coding_percentages.reindex(coding_types, fill_value=0)

    isoform_counts = isoform_counts.reindex(coding_types, fill_value=0)
    isoform_percentages = isoform_percentages.reindex(coding_types, fill_value=0)

    coding_counts_df = pd.DataFrame(coding_counts).reset_index()
    coding_counts_df['gene_name'] = gene_name

    coding_percentages_df = pd.DataFrame(coding_percentages).reset_index()
    coding_percentages_df['gene_name'] = gene_name

    isoform_counts_df = pd.DataFrame(isoform_counts).reset_index()
    isoform_counts_df['gene_name'] = gene_name

    isoform_percentages_df = pd.DataFrame(isoform_percentages).reset_index()
    isoform_percentages_df['gene_name'] = gene_name

    df_fl_assoc_counts_by_coding = pd.concat([df_fl_assoc_counts_by_coding, coding_counts_df])
    df_fl_assoc_percentages_by_coding = pd.concat([df_fl_assoc_percentages_by_coding, coding_percentages_df])
    df_isoform_counts_by_coding = pd.concat([df_isoform_counts_by_coding, isoform_counts_df])
    df_isoform_percentages_by_coding = pd.concat([df_isoform_percentages_by_coding, isoform_percentages_df])

filename = "R1_"
# Save the summary DataFrames to separate CSV files
df_fl_assoc_counts_by_level.to_csv(filename + 'fl_assoc_counts_by_level.csv', index=False)
df_fl_assoc_percentages_by_level.to_csv(filename + 'fl_assoc_percentages_by_level.csv', index=False)
df_isoform_counts_by_level.to_csv(filename + 'isoform_counts_by_level.csv', index=False)
df_isoform_percentages_by_level.to_csv(filename + 'isoform_percentages_by_level.csv', index=False)

df_fl_assoc_counts_by_coding.to_csv(filename + 'fl_assoc_counts_by_coding.csv', index=False)
df_fl_assoc_percentages_by_coding.to_csv(filename + 'fl_assoc_percentages_by_coding.csv', index=False)
df_isoform_counts_by_coding.to_csv(filename + 'isoform_counts_by_coding.csv', index=False)
df_isoform_percentages_by_coding.to_csv(filename + 'isoform_percentages_by_coding.csv', index=False)
