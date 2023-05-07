# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Classify NIS genes into functional subtypes and predict structural features
"""
from typing import Iterable, Dict, List

from antismash.common import fasta, path, subprocessing
from antismash.common.hmmscan_refinement import HMMResult, refine_hmmscan_results
from antismash.common.secmet import Record, CDSFeature
from antismash.common.utils import get_hmm_lengths
from antismash.modules.siderophores.results import SiderophoreResults


def run_nis_hmmscan(cds_features: Iterable[CDSFeature]) -> Dict[str, List[HMMResult]]:
    """ Runs hmmscan for type II PKS proteins on the given CDSFeatures

        Arguments:
            cds_features: CDSs that should be hmmscanned

        Returns:
            a dictionary of key: cds and value: list of HMMResults, for hmmscan results of the cluster
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "nis.hmm")
    # TODO: consider using TCs, where an NIS might have no hit
    hmm_results = subprocessing.run_hmmscan(hmm_file, cluster_fasta, opts=['-T 20'])  # opts=['--cut_tc']
    hmm_lengths = get_hmm_lengths(hmm_file)
    return refine_hmmscan_results(hmm_results, hmm_lengths)



def run_siderophore_analysis(record: Record) -> SiderophoreResults:
    """ Runs the siderophore analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated SiderophoreResults object
    """
    results = SiderophoreResults(record.id)
    for cluster in record.get_protoclusters():
        if cluster.product != "NI-siderophore":
            continue

        hmm_results_by_name = run_nis_hmmscan(cluster.definition_cdses)
        hmm_results_by_cds = {record.get_cds_by_name(name): hits for name, hits in hmm_results_by_name.items()}
