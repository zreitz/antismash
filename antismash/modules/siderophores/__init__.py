# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Based on functional subtypes of NI-siderophore synthetases, predict
    structural features of the produced siderophores
"""

# start with standard library imports
import logging
from typing import Any, Dict, List, Optional

# then any imports from external modules, e.g. biopython, if relevant

# then any imports from antismash
from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.common import hmmer

from .siderophore_analysis import run_siderophore_analysis


# then any local file imports, e.g. from .somefile import..., if relevant

NAME = "siderophores"
SHORT_DESCRIPTION = "Predict NI-siderophore structural features"


# define a results class, this is important as adding information to the record
# during analysis will cause issues
# for detailed examples, see any of the other analysis modules' implementations
class SiderophoreAnalysisResults(ModuleResults):
    """ Example results class for the analysis module template """  # TODO
    schema_version = 1  # when the data format in the results changes, this needs to be incremented

    # define whatever construction arguments are needed, record_id is required by the superclass
    # it's good to keep any command line option values here to know when they're changed for --reuse-results
    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.some_other_information = []  # this could be added to during analysis

    # implement a conversion to a JSON-compatible dictionary
    # all elements must one of: str, int, float, list, or a dict of those types (this can recurse)
    def to_json(self) -> Dict[str, Any]:
        """ Constructs a JSON representation of this instance """

        return {
            "schema_version": self.schema_version,
            "other": [str(item) for item in self.some_other_information],  # an example only
        }

    # once _all_ analysis modules have completed, their results are added with this method
    # adding to the record during the analysis will cause issues
    def add_to_record(self, record: Record) -> None:
        """ Adds the analysis results to the record """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
        # any results would be added here  # TODO
        # for an example of new features, see antismash.modules.tta
        # for an example of qualifiers, see antismash.modules.t2pks
        # any new feature types or qualifiers would be implemented in antismash.common.secmet,
        #   and would need to be able to be converted to and from biopython's SeqFeature without loss
        raise NotImplementedError()  # remove this when completed

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
        for other in json["other"]:
            results.some_other_information.append(other)

        return results


def get_arguments() -> ModuleArgs:
    """ Builds any commandline argument constructs that may be required

        Returns:
            an empty or populated ModuleArgs instance
    """
    # construct the argument group, with section and prefix
    # the prefix will be enforced for all command line options for the module
    args = ModuleArgs('Additional analysis', 'siderophore')

    # an example toggle to turn on your analysis, if not set to always be enabled
    # can also be used to turn on/off extra features of your analysis
    args.add_analysis_toggle('--siderophore',     # the option as it appears on the command line
                             dest='siderophore',  # the storage location in the antismash Config object
                             default=False,             # disabled by default
                             action='store_true',       # enabled if --siderophore-analysis is given on the commandline
                             # and finally, text to show when the user runs with --help
                             help="Predict NI-siderophore structural features.")

    return args


def check_options(options: ConfigType) -> List[str]:
    """ Checks that the provided options are compatible with each other

        Arguments:
            options: the current antismash config object

        Returns:
            a list of strings, each string being an issue with the given options
    """
    return []


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check the prerequisites.
            hmmscan: domain detection
            HMMs: ni-siderophore.hmm

        Returns:
            a list of strings describing any errors, if they occurred
    """
    failure_messages = []
    for binary_name in ['hmmscan', "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable: {binary_name}")

    failure_messages.extend(prepare_data(logging_only=True))

    return failure_messages


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    database = path.get_full_path(__file__, "data", "nis.hmm")
    return hmmer.ensure_database_pressed(database, return_not_raise=logging_only)


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled by the given options """
    return options.siderophore or not options.minimal


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[SiderophoreAnalysisResults]:
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    # if there isn't anything to work with, just return None
    if not previous:
        return None
    return SiderophoreAnalysisResults.from_json(previous, record)


def run_on_record(record: Record, results: SiderophoreAnalysisResults, _options: ConfigType) -> SiderophoreAnalysisResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation
    """
    # after a safety check that the results are the correct ones for the record, return them
    if isinstance(results, SiderophoreAnalysisResults) and results.record_id == record.id:
        return results
    # otherwise run the actual analysis and generate a results instance with your analysis results
    results = run_siderophore_analysis(record)
    # and return it
    return results
