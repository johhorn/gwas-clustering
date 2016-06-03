import sys

"""
Takes input file of SNPs and their LD variants on the following format:
<chrnum>    <ldPos> <ldSNP> <tagPos>    <tagSNP>    <r2>

Creates a number of linked point tracks, keeping variants with ld scores above
the given r2
"""


def add_edge(edge_dict, rsid, start, seqid, edge_rsid, r2):
    if rsid in edge_dict:
        if edge_rsid not in edge_dict[rsid]['edges']:
            edge_dict[rsid]['edges'][edge_rsid] = r2
    else:
        edge_dict[rsid] = {
            'start': start,
            'seqid': seqid,
            'edges': {edge_rsid: r2}
        }


def find_edges(SNPfile, edge_dict):

    for line in SNPfile:

        snp = line.strip().split('\t')

        tagSNP = snp[4]
        ldSNP = snp[2]

        if not tagSNP.startswith('rs') and not ldSNP.startswith('rs'):
            continue
        if tagSNP == ldSNP:
            continue

        chrnum = snp[0]
        tagPos = snp[3]
        ldPos = snp[1]
        r2 = float(snp[5])

        add_edge(edge_dict, rsid=tagSNP, start=tagPos, seqid=chrnum, edge_rsid=ldSNP, r2=r2)
        add_edge(edge_dict, rsid=ldSNP, start=ldPos, seqid=chrnum, edge_rsid=tagSNP, r2=r2)


def write_header_info(filename):
    filename.write('##gtrack version: 1.0\n')
    filename.write('##track type: linked points\n')
    filename.write('##undirected edges: True\n')
    filename.write('##edge weights: True\n')
    filename.write('##edge weight type: number\n')
    filename.write('##edge weight dimension: scalar\n')
    filename.write('##uninterrupted data lines: true\n')
    filename.write('##no overlapping elements: true\n')
    filename.write('##1-indexed: false\n')
    filename.write('###seqid\tstart\tid\tedges\n')

def create_linked_point_track(edge_dict, r2_filter):
    gtrack_file = open('linked_point_track_' + str(int(r2_filter*10)) + '.gtrack', 'w')
    write_header_info(gtrack_file)

    for rsid, snpInfo in edge_dict.items():
        start = snpInfo['start']
        seqid = 'chr' + str(snpInfo['seqid'])
        edges = snpInfo['edges']

        text_edges = ""
        separator = ""
        for snp in edges:
            if edges[snp] < r2_filter:
                continue
            else:
                edge = separator + snp + '=' + str(edges[snp])
                text_edges += edge
                separator = ';'

        if text_edges != "":
            gtrack_file.write(seqid + '\t' + start + '\t' + rsid + '\t' + text_edges + '\n')
    gtrack_file.close()


def main():
    try:
        snpFilename = sys.argv[1]
        r2 = float(sys.argv[2])
    except:
        print 'python', sys.argv[0], '<expanded snps file> <r2>'
        exit()

    SNPfile = open(snpFilename, 'r')

    edge_dict = {}
    find_edges(SNPfile, edge_dict)
    create_linked_point_track(edge_dict, r2)

    SNPfile.close()


if __name__ == '__main__':
    main()
