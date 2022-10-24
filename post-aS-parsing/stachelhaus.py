import csv


with open("nrpspredictor2_sids.tsv", 'r') as codef:
    codes = list(csv.reader(codef, delimiter = "\t"))

lookup = {}
for line in codes:
    code = line[1]
    aa = line[0].split("_")[0]
    if code not in lookup:
        lookup[code] = aa

with open("triserine.tsv", 'r') as resultf:
    results = list(csv.reader(resultf, delimiter = "\t"))

matches = []
for nrps in results[1:]:
    code = nrps[5].upper()
    bestaa = set()
    bestmatch = 0
    for known_code, aa in lookup.items():
        match = sum(known_code_ == code_ for known_code_, code_ in zip(known_code, code))
        if match > bestmatch:
            bestmatch = match
            bestaa = {aa,}
        elif match == bestmatch:
            bestaa.add(aa)
    nrps.extend(["/".join(bestaa), bestmatch])
    matches.append(nrps)

with open("triserine-re-pred.tsv", 'w') as outf:
    writer = csv.writer(outf, delimiter = "\t")
    for row in matches:
        writer.writerow(row)