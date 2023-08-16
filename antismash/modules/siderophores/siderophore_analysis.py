# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Classify NIS genes into functional subtypes and predict structural features
"""
from collections import defaultdict

from antismash.common import path, hmmer
from antismash.common.hmmer import HmmerResults
from antismash.common.secmet import Record
from antismash.modules.siderophores import SiderophoreAnalysisResults
from antismash.modules.siderophores.siderophore_classes import SiderophorePrediction



def identify_pathways(hmmer_results: HmmerResults):
    # Reshape to HMM focused
    hmms = {}
    for hit in hmmer_results.hits:
        hmms.setdefault(hit.domain, []).append(hit)

    # Certain HMMs are/aren't expected to co-occur. Warn the user otherwise.
    caveats = []
    # Caveat if there's an IucA_IucC hit that doesn't have a specific hit

    # Precursor biosynthesis

    # Core biosynthesis




def run_siderophore_analysis(record: Record) -> SiderophoreAnalysisResults:
    """ Runs the siderophore analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated SiderophoreAnalysisResults object
    """

    nis_db = path.get_full_path(__file__, "data", "nis.hmm")
    hmmer_results_list = []
    predictions = []
    for cluster in record.get_protoclusters():
        if cluster.product != "NI-siderophore":
            continue

        # Scan for biosynthetic (sub-)families
        hmmer_results = hmmer.run_hmmer(record,
                                        cluster.cds_children,
                                        max_evalue=0.1, min_score=10,
                                        database=nis_db,
                                        tool="siderophores")
        hmmer_results_list.append(hmmer_results)

        # Construct biosynthetic schemes
        predictions.append(SiderophorePrediction(cluster.get_protocluster_number(),
                                                 hmmer_results))

    results = SiderophoreAnalysisResults.from_hmmer_results(hmmer_results)
    return results
