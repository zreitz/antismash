## Given a directory of antismash results, perform analyses

import os
from sys import argv
import pandas as pd
import json

from tqdm import tqdm


def parse_results(asdir):
    record_infos = []
    for genome in tqdm(os.listdir(asdir)):
        try: # Lazy
            if not genome.startswith("GCF"):
                continue

            genome_dir = os.path.join(asdir, genome)
            regions = (r for r in os.listdir(genome_dir) if "region" in r and not r.startswith("GCF"))
            if not regions:
                continue

            jsonpath = os.path.join(genome_dir, [s for s in os.listdir(genome_dir) if ".json" in s][0])
            with open(jsonpath, 'r') as jsonf:
                data = json.load(jsonf)

            # Check for 2,3-DHB
            has_cat = False
            for record in data['records']:
                for area in record.get("areas"):
                    if "23-DHB" in area.get("products"):
                        has_cat = True
                        break
                else:
                    continue
                break
            if not has_cat:
                continue

            for record in (r for r in data['records'] if len(r['areas']) > 0):
                cds_results = record['modules']['antismash.detection.nrps_pks_domains']['cds_results']
                for cds, results in cds_results.items():
                    domains = [h['hit_id'] for h in results['domain_hmms']]
                    if not len(domains) in (7, 8, 9): # With/without E and TIGR01720
                        continue
                    if not (domains[:3] == ['Condensation_Starter', 'AMP-binding', 'PCP']
                            and domains[-3:] == ['AMP-binding', 'PCP', 'Thioesterase']):
                        continue
                    domain_preds = record['modules']['antismash.modules.nrps_pks']['domain_predictions']
                    this_cds_As = [p for p in domain_preds.items() if cds in p[0] and "AMP-binding" in p[0]]
                    if len(this_cds_As) != 2:
                        continue
                    first = [i[1] for i in this_cds_As if i[0].endswith("1")][0]
                    second = [i[1] for i in this_cds_As if i[0].endswith("2")][0]
                    cds_dict = {"Assembly": genome,
                                "Strain": record['annotations']['organism'],
                                "NRPS": cds,
                                "A1.Prediction": "/".join(first['NRPSPredictor2']['stachelhaus_predictions']),
                                "A1.Stachelhaus": first['NRPSPredictor2']['stachelhaus_seq'],
                                "A2.Prediction": "/".join(second['NRPSPredictor2']['stachelhaus_predictions']),
                                "A2.Stachelhaus": second['NRPSPredictor2']['stachelhaus_seq'],
                                "Edomain": "Epimerization" in domains}
                    if cds_dict["A2.Prediction"] != "ser":
                        continue
                    record_infos.append(cds_dict)
        except IndexError:
            print(genome)
            record_infos.append({"Assembly": genome, "NRPS": "ERROR"})



    # Convert to dataframe (except last column)
    df = pd.DataFrame((r.values() for r in record_infos),
                      columns = ("Assembly", "Strain", "NRPS", "A1.Prediction", "A1.Stachelhaus", "A2.Prediction", "A2.Stachelhaus", "Edomain"))


    outpath = os.path.join(asdir, "triserine.tsv")
    print("Writing score table to", outpath)
    with open(outpath, 'w') as outf:
        df.to_csv(outf, sep = "\t")


if __name__ == "__main__":
    asdir = argv[1]
    parse_results(asdir)
