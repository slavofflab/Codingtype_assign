import pandas as pd
import re
from collections import defaultdict
from pyfaidx import Fasta
from Bio.Seq import Seq

# === 1. 加载 iORF 数据 ===
iorf_df = pd.read_csv("V35_iORF_NCBI.csv")
#iorf_df = pd.read_csv("V35_iORF - Copy.csv")
iorf_df.columns = [c.strip().lower() for c in iorf_df.columns]
iorf_df = iorf_df.rename(columns={"chrm": "chr"})
iorf_dict = dict(zip(iorf_df["orf_name"], iorf_df["gene_name"]))

# === 2. 加载 GTF 数据 ===
gtf_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_csv("E:/DataBase/hg38_v38/GRCh38_latest_genomic.gff", sep="\t", comment="#", header=None, names=gtf_cols, dtype=str)
#gtf = pd.read_csv("test - Copy.gtf", sep="\t", comment="#", header=None, names=gtf_cols, dtype=str)
genome = Fasta("E:/DataBase/pri_hg38.fa")

gtf = gtf[gtf["seqname"].str.startswith("NC_")].copy()
# 提取 key=value 的字段，去除版本号
def extract_attr(attr, key):
    match = re.search(rf'{key}=([^;]+)', attr)
    return match.group(1).split(".")[0] if match else None

# 提取 transcript_id 和 gene_id 所需字段
gtf["ID"] = gtf["attribute"].apply(lambda x: extract_attr(x, "ID"))
gtf["Parent"] = gtf["attribute"].apply(lambda x: extract_attr(x, "Parent"))
gtf["start"] = gtf["start"].astype(int)
gtf["end"] = gtf["end"].astype(int)

# 构建 transcript_id ↔ gene_id 映射（只使用 mRNA / transcript 行）
transcript_to_gene = {}
for _, row in gtf.iterrows():
    if row["feature"] in ["mRNA", "transcript"]:
        tid = extract_attr(row["attribute"], "ID")
        gid = extract_attr(row["attribute"], "gene")
        if tid and gid:
            tid = tid.split(".")[0]
            gid = gid.split(".")[0]
            transcript_to_gene[tid] = gid

# 设置统一 transcript_id 字段（mRNA/transcript 用 ID，其余用 Parent）
gtf["transcript_id"] = gtf.apply(
    lambda row: extract_attr(row["attribute"], "ID") if row["feature"] in ["mRNA", "transcript"]
    else extract_attr(row["attribute"], "Parent"),
    axis=1
)
gtf["transcript_id"] = gtf["transcript_id"].str.split(".").str[0]  # 去版本号


# === 3. 构建 transcript 特征：CDS / start_codon / exon / strand ===
transcript_features = defaultdict(lambda: {"CDS": [], "exon": [], "strand": None})
gene_to_transcripts = defaultdict(list)

for _, row in gtf.iterrows():
    tid = row["transcript_id"]
    if pd.isna(tid):
        continue

    gid = transcript_to_gene.get(tid)
    if gid is None:
        continue

    # gene_to_transcripts: 只在 mRNA / transcript 行记录一次
    if row["feature"] in ["mRNA", "transcript"]:
        if tid not in gene_to_transcripts[gid]:
            gene_to_transcripts[gid].append(tid)

    # transcript_features
    if row["feature"] == "CDS":
        transcript_features[tid]["CDS"].append((row["start"], row["end"]))
    elif row["feature"] == "exon":
        transcript_features[tid]["exon"].append((row["start"], row["end"]))
    transcript_features[tid]["strand"] = row["strand"]

# === 4. 工具函数 ===
def to_pos_list(intervals, strand):
    pos = []
    for start, end in intervals:
        pos.extend(range(start, end + 1))  # always from low to high

    # 如果是负链，需要反向整个 pos 列表
    if strand == "-":
        pos.reverse()
    
    return pos


def get_iorf_intervals(row):
    starts = list(map(int, str(row["starts"]).split(";")))
    ends = list(map(int, str(row["ends"]).split(";")))
    return list(zip(starts, ends))

def is_sublist(small, big):
    #Check if `small` list appears as a consecutive sublist in `big` list.
    n = len(small)
    for i in range(len(big) - n + 1):
        if big[i:i + n] == small:
            return True
    return False

def merge_codon_parts_simple(codon_intervals):
    if len(codon_intervals) <= 1:
        return codon_intervals
    else:
        starts = [s for s, e in codon_intervals]
        ends = [e for s, e in codon_intervals]
        return [(min(starts), max(ends))]
def codon_sort_key(codon, strand):
    return codon[0] if strand == "+" else -codon[1]
def normalize_codon_orientation(codon, strand):
    start, end = codon
    return (start, end) if strand == "+" else (end, start)
def coframe_stats(iorf_pos_list: list[int], cds_pos_list: list[int]) -> tuple[int, float]:
    cds_idx_by_pos = {coord: idx for idx, coord in enumerate(cds_pos_list)}

    coframe_count = 0
    for iorf_idx, coord in enumerate(iorf_pos_list):
        cds_idx = cds_idx_by_pos.get(coord)
        if cds_idx is not None and (iorf_idx % 3) == (cds_idx % 3):
            coframe_count += 1

    iorf_len = len(iorf_pos_list)
    coframe_ratio_all = (coframe_count / iorf_len) if iorf_len else 0.0

    return coframe_count, coframe_ratio_all


# === 5. 分析每个 iORF ===
results = []
for _, row in iorf_df.iterrows():
    gid = row["gene_name"]
    iorf_strand = row["strand"]
    iorf_id = row["orf_name"]

    if pd.isna(gid) or gid not in gene_to_transcripts:
        print(gid)
        continue

    iorf_intervals = get_iorf_intervals(row)
    iorf_intervals = sorted(iorf_intervals, key=lambda x: x[0])
    iorf_pos_list = to_pos_list(iorf_intervals, iorf_strand)
    
    # === 提取该基因所有转录本的 start_codon / stop_codon 区段（通过 CDS 起止推断） ===
    all_start_codons = []
    all_stop_codons = []

    for tid in gene_to_transcripts[gid]:
        t_feats = transcript_features.get(tid)
        if not t_feats or not t_feats["CDS"]:
            continue

        strand = t_feats["strand"]
        cds_sorted = sorted(t_feats["CDS"], key=lambda x: x[0])

        if strand == "+":
            start_codon = (cds_sorted[0][0], cds_sorted[0][0] + 2)
            stop_codon = (cds_sorted[-1][1] - 2, cds_sorted[-1][1])
        else:
            start_codon = (cds_sorted[-1][1] - 2, cds_sorted[-1][1])
            stop_codon = (cds_sorted[0][0], cds_sorted[0][0] + 2)

        # 写入 transcript_features，方便后续使用
        t_feats["start_codon"] = [start_codon]
        t_feats["stop_codon"] = [stop_codon]

        all_start_codons.append(start_codon)
        all_stop_codons.append(stop_codon)

    # 去重 + 按照 strand 排序（统一方向为 5′→3′，即 iORF 所在链的方向）
    all_start_codons = sorted(set(all_start_codons), key=lambda x: x[0], reverse=(iorf_strand == "-"))
    all_stop_codons = sorted(set(all_stop_codons), key=lambda x: x[0], reverse=(iorf_strand == "-"))
    all_start_codons = [normalize_codon_orientation(c, iorf_strand) for c in all_start_codons]
    all_stop_codons  = [normalize_codon_orientation(c, iorf_strand) for c in all_stop_codons]

    # === 遍历该基因的每个 transcript 进行后续分析 ===
    for tid in set(gene_to_transcripts[gid]):
        t_feats = transcript_features.get(tid)
        if not t_feats:
            continue

        cds_intervals = sorted(t_feats["CDS"], key=lambda x: x[0])
        exon_intervals = sorted(t_feats["exon"], key=lambda x: x[0])
        start_codons = t_feats.get("start_codon", [])
        stop_codon_sites = t_feats.get("stop_codon", [])
        cds_strand = t_feats["strand"]
        chrom = "chr" + str(row["chr"])  # 统一格式


        exon_pos_list = to_pos_list(exon_intervals, cds_strand) if exon_intervals else []
        cds_pos_list = to_pos_list(cds_intervals, cds_strand) if cds_intervals else []

        # === 提取拼接转录本序列 ===
        try:
            exon_seqs = [genome[chrom][start - 1:end].seq for start, end in exon_intervals]
            transcript_seq = "".join(exon_seqs)
            if cds_strand == "-":
                transcript_seq = str(Seq(transcript_seq).reverse_complement())
        except Exception as e:
            print(f"[ERROR] {tid} {chrom} sequence failed: {e}")
            transcript_seq = None

        has_start = "yes" if start_codons else "no"
        has_stop_annotation = bool(stop_codon_sites)
        has_cds = bool(cds_intervals) and bool(start_codons)
        transcript_length = len(exon_pos_list)
        iorf_length = len(iorf_pos_list)

        # iORF 是否为 transcript 的子区段
        iorf_in_transcript = "yes" if is_sublist(iorf_pos_list, exon_pos_list) else "no"

        # iORF 最后21nt 是否为 transcript 的子区段
        last_n = 21
        last21nt_of_iorf_in_transcript = "yes" if is_sublist(iorf_pos_list[-last_n:], exon_pos_list) else "no"

        if iorf_in_transcript == "yes":
            iorf_start_in_transcript = exon_pos_list.index(iorf_pos_list[0]) + 1  # 1-based
            iorf_end_in_transcript = exon_pos_list.index(iorf_pos_list[-1]) + 1
        else:
            iorf_start_in_transcript = iorf_end_in_transcript = None

        # Upsteam annotated start codon
        has_up_start_annotation = False

        if iorf_in_transcript == "yes":
            for start, _ in all_start_codons:
                try:
                    codon_start_tx_pos = exon_pos_list.index(start) + 1  # 1-based transcript coordinate
                    if codon_start_tx_pos < iorf_start_in_transcript:
                        has_up_start_annotation = True
                        break  # 找到一个即可
                except ValueError:
                    continue  # start codon 不在这个 transcript 上

        if has_cds:
            cds_start_in_transcript = exon_pos_list.index(cds_pos_list[0]) + 1
            cds_end_in_transcript = exon_pos_list.index(cds_pos_list[-1]) + 1
            # === 尝试提取 cds_end_in_transcript 后面三个碱基并判断是否为 stop codon ===
            try:
                # 使用 0-based 索引定位 CDS 结束位置
                cds_end_index_0_based = exon_pos_list.index(cds_pos_list[-1])
                downstream_codon = transcript_seq[cds_end_index_0_based:cds_end_index_0_based+3]
                if downstream_codon in {"TAA", "TAG", "TGA"}:
                    has_stop_annotation = True
            except Exception:
                pass
        else:
            cds_start_in_transcript = cds_end_in_transcript = None

        # 计算距离 iORF 结束位点与最后一个 exon 起始位点在 transcript 上的距离
        if iorf_in_transcript == "yes":
            if cds_strand == "+":
                last_exon_start_genomic = exon_intervals[-1][0]
                last_exon_start_in_transcript = exon_pos_list.index(last_exon_start_genomic)
            else:
                last_exon_start_genomic = exon_intervals[0][1]
                last_exon_start_in_transcript = exon_pos_list.index(last_exon_start_genomic)
            iorf_end_idx = exon_pos_list.index(iorf_pos_list[-1])
            distance_to_last_exon_start = last_exon_start_in_transcript - iorf_end_idx
        else:
            distance_to_last_exon_start = None

        # === 提取 CDS → iORF 的序列并判断是否存在 in-frame stop codon ===
        early_stop = "NA"
        cds_to_iorf_seq = ""
        upstream_inframe_atg = "NA"

        if transcript_seq and has_cds and iorf_in_transcript == "yes":
            try:
                cds_start_idx = exon_pos_list.index(cds_pos_list[0])  # 0-based
                iorf_start_idx = exon_pos_list.index(iorf_pos_list[0])  # 0-based

                # === New: upstream scan in the iORF reading frame (from CDS start toward 5') ===
                seq_upper = transcript_seq.upper()
                stop_codon_set = {"TAA", "TAG", "TGA"}

                # align to iORF frame at or before cds_start_idx
                offset = (cds_start_idx - iorf_start_idx) % 3
                start_pos = cds_start_idx - offset  # nearest ≤ cds_start_idx in iORF frame

                upstream_inframe_atg = "no"
                upstream_atg_pos = None
                codons_to_iorf = None
                codons_to_cds  = None

                # exclude the codon at cds_start_idx itself; scan strictly upstream in iORF frame
                for pos in range(start_pos - 3, -1, -3):
                    codon = seq_upper[pos:pos+3]
                    if len(codon) < 3:
                        break
                    if codon in stop_codon_set:
                        # stop encountered before any in-frame ATG → terminate upstream search
                        break
                    if codon == "ATG":
                        upstream_inframe_atg = "yes"
                        upstream_atg_pos = pos
                        codons_to_iorf = (iorf_start_idx - pos) // 3
                        codons_to_cds  = (cds_start_idx - pos) // 3
                        break
                # === End of new part ===

                if iorf_start_idx > cds_start_idx:
                    cds_to_iorf_seq = transcript_seq[cds_start_idx:iorf_start_idx]
                    frame_offset = (iorf_start_idx - cds_start_idx) % 3
                    trimmed_seq = cds_to_iorf_seq[frame_offset:]
                    # keep your original early_stop logic unchanged
                    has_early_stop_codon = any(
                        trimmed_seq[i:i+3] in stop_codon_set
                        for i in range(0, len(trimmed_seq) - 2, 3)
                    )
                    early_stop = "yes" if has_early_stop_codon else "no"
                else:
                    early_stop = "not_after_CDS"

            except Exception as e:
                print(f"[ERROR parsing iORF→CDS seq: {tid}, {iorf_id}] {e}")

        # === 分类判断 ===
        if iorf_in_transcript == "no":
            rel_pos = "not_in_transcript"
        elif not has_cds and has_up_start_annotation:
            rel_pos = "iORF in no stop"
        elif not has_cds:
            rel_pos = "ncRNA_ORF"
        else:
            if iorf_end_in_transcript < cds_start_in_transcript:
                rel_pos = "uORF"
            elif (
                iorf_start_in_transcript < cds_start_in_transcript
                and iorf_end_in_transcript >= cds_start_in_transcript
                and iorf_end_in_transcript < cds_end_in_transcript
            ):
                rel_pos = "uoORF"
            elif (
                iorf_start_in_transcript < cds_start_in_transcript
                and iorf_end_in_transcript >= cds_end_in_transcript
            ):
                rel_pos = "overprintORF"
            elif (
                iorf_start_in_transcript == cds_start_in_transcript
                # and iorf_end_in_transcript == cds_end_in_transcript 
                # → This was ignored because new ORFs sometimes lack stop codon annotation
            ):
                rel_pos = "CDS"
            elif (
                iorf_start_in_transcript > cds_start_in_transcript
                and iorf_end_in_transcript <= cds_end_in_transcript
            ):
                rel_pos = "iORF"
            elif (
                iorf_start_in_transcript > cds_start_in_transcript
                and (
                    (iorf_end_in_transcript - 3 == cds_end_in_transcript) or
                    (iorf_end_in_transcript == cds_end_in_transcript)
                )
            ):
                rel_pos = "C-terminal isoform"
            elif (
                iorf_start_in_transcript > cds_start_in_transcript
                and iorf_end_in_transcript > cds_end_in_transcript
            ):
                rel_pos = "doORF"
            elif iorf_start_in_transcript > cds_end_in_transcript:
                rel_pos = "dORF"
            else:
                rel_pos = "other"

        # === 判断是否与任何 stop codon 注释共用 C 端 ===
        if last21nt_of_iorf_in_transcript == "yes":
            try:
                iorf_end_in_transcript = exon_pos_list.index(iorf_pos_list[-1]) + 1  # 1-based transcript coordinate
            except ValueError:
                same_c_term = "NA"  # 如果 iORF end 不在 exon_pos_list 中
            else:
                same_c_term = "no"
                for stop_start, stop_end in all_stop_codons:
                    try:
                        stop_end_in_tx = exon_pos_list.index(stop_end) + 1
                        if iorf_end_in_transcript == stop_end_in_tx or iorf_end_in_transcript - 3 == stop_end_in_tx:
                            same_c_term = "yes"
                            break
                    except ValueError:
                        continue  # stop codon 不在该 transcript 上
        else:
            same_c_term = "NA"

        #coframe count
        co_count, co_ratio = "NA", "NA"
        if has_cds:
            co_count, co_ratio = coframe_stats(iorf_pos_list, cds_pos_list)

        results.append({
            "iorf_id": iorf_id,
            "gene_name": iorf_dict.get(iorf_id, None),
            "gene_id": gid,
            "transcript_id": tid.replace("rna-", ""),
            "iorf_strand": iorf_strand,
            "cds_strand": cds_strand,
            "iorf_in_transcript": iorf_in_transcript,
            "has_upstream_start_annotation": has_up_start_annotation,
            "has_start_codon": has_start,
            "ORF_biotype": rel_pos,
            "transcript_length": transcript_length,
            "iorf_length": iorf_length,
            "iorf_start_in_transcript": iorf_start_in_transcript,
            "iorf_end_in_transcript": iorf_end_in_transcript,
            "cds_start_in_transcript": cds_start_in_transcript,
            "cds_end_in_transcript": cds_end_in_transcript,
            "distance_to_last_exon_start": distance_to_last_exon_start,
            "iORF_CDS_same_C_term": same_c_term,
            "CDS_coframe_count":co_count,
            "CDS_coframe_ratio":co_ratio,
            "early_stop_in_frame": early_stop,
            "upstream_inframe_atg": upstream_inframe_atg,
            "cds_to_iorf_seq": cds_to_iorf_seq,
            "transcript_seq": transcript_seq
        })

# === 6. 保存结果 ===
result_df = pd.DataFrame(results)
result_df.to_csv("iORF_transcript_analysis_NCBI.csv", index=False)
print("✅ 分析完成，结果已保存至 iORF_transcript_analysis.csv")

