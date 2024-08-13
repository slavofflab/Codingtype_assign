def attr_interpreter(attr):
    attr_list = {}
    # Split the string and remove extra spaces
    attr_s = attr.split()
    l = int(len(attr_s) / 2)
    for i in range(l):
        key = attr_s[2 * i]
        value = attr_s[2 * i + 1].strip(";").strip("\"")
        attr_list[key] = value
    return attr_list

class gtf_line:
    def __init__(self, g_line):
        line_tem = g_line.strip().split("\t")
        if len(line_tem) != 9:
            raise ValueError("GTF line does not have 9 fields.")
        self.chrom = line_tem[0]
        self.source = line_tem[1]
        self.type = line_tem[2]
        self.start = int(line_tem[3])
        self.end = int(line_tem[4])
        self.score = line_tem[5]
        self.strand = line_tem[6]
        self.phase = line_tem[7]
        self.attributes = attr_interpreter(line_tem[8])
    def __str__(self):
        return f"GTF Line: {self.chrom}, {self.source}, {self.type}, {self.start}, {self.end}, {self.strand}, {self.phase}, {self.attributes}"
