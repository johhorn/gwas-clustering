import requests, sys, time

inFile = sys.argv[1]
r2_threshold = sys.argv[2]

rsids = open(inFile)
snps = []
for line in rsids:
	rsid = line.strip()
	if rsid.startswith('rs'):
		snps.append(rsid)

def sortRsidTuple(rsid1, rsid2):
    id1 = int(rsid1[2:])
    id2 = int(rsid2[2:])
    return (rsid1, rsid2) if id1 < id2 else (rsid2, rsid1)


def addEdge(r2graph, rsid1, rsid2, r2):
        key = sortRsidTuple(rsid1, rsid2)
        if key not in r2graph:
            r2graph[key] = r2

server = "https://rest.ensembl.org"
api = "/ld/human/"
population = "?population_name=1000GENOMES:phase_3:CEU"
r2 = ";r2=" + str(r2_threshold)
windows= ";window_size=500"
r2graph = {}

for snp in snps:
	print 'extracting for snp', snp

	ext = api + snp + population + r2 + windows
	r = requests.get(server+ext, headers={"Content-Type" : "application/json"})
	decoded = r.json()

	if r.status_code == 429:
		print 'rate reached'
		seconds = r.headers('Retry-After')
		print seconds
		time.sleep(seconds + 10)
	elif r.status_code == 400:
		print 'Bad request'
		if 'error' in decoded:
			print decoded['error']
	elif not r.ok:
  		r.raise_for_status()
  		sys.exit()
  	else:
		for response in decoded:
			r2 = response['r2']
			v1 = response['variation1']
			v2 = response['variation2']
			addEdge(r2graph, v1, v2, r2)

print r2graph