# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classes for the siderophore analysis module. """

from typing import Any, Dict, List, Optional, TypeVar

from Bio.SeqFeature import FeatureLocation

from antismash.common.hmmer import HmmerResults, HmmerHit
from antismash.common.secmet import Record, AntismashDomain
from antismash.common.secmet.locations import Location, location_from_string

T = TypeVar("T", bound="SiderophoreDomain")
TOOL = "siderophores"

DOMAIN_TYPE_MAPPING = {
    "TypeA_Lys": "A",
    "TypeA_SerDap": "A",
    "TypeA_diamine": "A",
    "TypeA_polyamine": "A",
    "TypeA_diamine_Strep": "A",
    "TypeA_SerDap_Strep": "A",

    "TypeAprime": "Aprime",
    "TypeAprime_fungal": "Aprime",

    "TypeB": "B",
    "TypeB_Strep": "B",

    "TypeC_DapDab": "C",
    "TypeC_Lys": "C",
    "TypeC_diamine": "C",
    "TypeC_polyamine": "C",
    "TypeC_diamine_Strep": "C",

    "TypeCprime": "Cprime",

    "Citrate_synth": "precursor",
    "Dab-synth": "precursor",
    "SbnA": "precursor",
    "SbnB": "precursor",

    "Pyridoxal_deC": "amine",
    "Lys_Orn_oxgnase": "amine",
    "Acetyltransf_8": "amine",
    "Penicil_amidase": "amine",
    "Acetyltransf_1": "amine",
    "Orn_Arg_deC_N": "amine",

    "DHS-dehydratase": "34dhb",
    "34-DHB-transferase": "34dhb",

    # Not used for predictions:
    # FhuF
}



class SiderophoreAnalysisResults(HmmerResults):
    """ Holds the results of siderophore analysis for a record """
    schema_version = 1

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    def __init__(self, record_id: str, evalue: float, score: float, hits: List[HmmerHit]) -> None:
        super().__init__(record_id, evalue, score, "nis.hmm", "Siderophores", hits)
        # self.some_other_information = []  # this could be added to during analysis


    # implement a conversion to a JSON-compatible dictionary
    # all elements must one of: str, int, float, list, or a dict of those types (this can recurse)
    # TODO: for now using HmmerResults.to_json()
    # def to_json(self) -> Dict[str, Any]:
    #     """ Constructs a JSON representation of this instance """
    #
    #     return {
    #         "schema_version": self.schema_version,
    #         "other": [str(item) for item in self.some_other_information],  # an example only
    #     }

    # once _all_ analysis modules have completed, their results are added with this method
    # adding to the record during the analysis will cause issues
    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")

        for i, hit in enumerate(self.hits):
            protein_location = FeatureLocation(hit.protein_start, hit.protein_end)
            nis_feature = AntismashDomain(location_from_string(hit.location),
                                          tool=TOOL,
                                          protein_location=protein_location,
                                          locus_tag=hit.locus_tag,
                                          domain=hit.domain)
            nis_feature.detection = "hmmscan"
            nis_feature.domain_id = f"{self.tool}_{nis_feature.locus_tag}_{i + 1:04d}"
            record.add_feature(nis_feature)

    # implement a conversion from the JSON-compatible data returned by to_json()
    # this allows --results-reuse to avoid running the module again if not neccessary
    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["SiderophoreAnalysisResults"]:
        """ Constructs a new results instance from a JSON format and the
            original record analysed. # TODO
        """
        # check that the previous data version is the same as current, if not, discard the results
        if json["schema_version"] != SiderophoreAnalysisResults.schema_version:
            return None

        # the exact reconstruction depends on what details are stored
        # to match the conversion to JSON that would be:
        results = SiderophoreAnalysisResults(json["record_id"])
        # for other in json["other"]:
        #     results.some_other_information.append(other)

        return results

    @classmethod
    def from_hmmer_results(cls, hmmer_results: HmmerResults) -> "SiderophoreResults":
        """ Convert HmmerResults into SiderophoreResults. """
        return cls(hmmer_results.record_id, hmmer_results.evalue, hmmer_results.score, hmmer_results.hits)



class SiderophorePrediction:
    schema_version = 1

    def __init__(self, protocluster_nr: int, hmmer_results: HmmerResults):
        self.protocluster_number = protocluster_nr
        self.precursors = {}
        self.core = []
        self.caveats = []
        self.identify_pathways(hmmer_results)

    def identify_pathways(self, hmmer_results: HmmerResults):
        # Reshape to HMM focused
        domains = {}
        for hit in hmmer_results.hits:
            #domain_type = DOMAIN_TYPE_MAPPING.get(hit.domain, "other")
            domains.setdefault(hit.domain, []).append(hit)

        # Precursors
        if "Citrate_synth" in domains:
            self.precursors["Citrate"] = domains["Citrate_synth"]
        if "Dab-synth" in domains:
            self.precursors["Citrate"] = domains["Dab-synth"]


class Reaction(List):
    def __init__(self, reactants: list[str],
                 enzymes: list[HmmerHit],
                 products: list[str],
                 underarrow: list[str] = None):
        super().__init__((reactants, [enzymes, underarrow], products))


    # By default
    def add_step(self, enzymes:list[HmmerHit], products: list[str],
                 underarrow: list[str] = None, left: bool = False):
        pass