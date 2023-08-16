# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Based on functional subtypes of NI-siderophore synthetases, predict
    structural features of the produced siderophores
"""

# start with standard library imports
import logging
from typing import Any, Dict, List, Optional, TypeVar

# then any imports from external modules, e.g. biopython, if relevant

# then any imports from antismash
from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.common import hmmer

# then any local file imports, e.g. from .somefile import..., if relevant
from .siderophore_classes import SiderophoreAnalysisResults
from .siderophore_analysis import run_siderophore_analysis
from ...common.secmet.features.antismash_domain import register_asdomain_variant

NAME = "siderophores"
SHORT_DESCRIPTION = "Predict NI-siderophore structural features"

# T = TypeVar("T", bound="SiderophoreDomain")
# TOOL = "siderophores"


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


#register_asdomain_variant(TOOL, TIGRDomain)