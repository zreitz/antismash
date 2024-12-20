# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generic sideloading """

import itertools
import logging
from typing import List, Optional

from antismash.common import path
from antismash.common.errors import AntismashInputError
from antismash.common.secmet import Record

from .data_structures import (
    ProtoclusterAnnotation,
    SideloadSimple,
    SideloadedResults,
    SubRegionAnnotation,
    Tool,
)
from .loader import load_validated_json

_SCHEMA_FILE = path.get_full_path(__file__, "schemas", "general", "schema.json")


def load_single_record_annotations(annotation_files: List[str], record: Record,
                                   manual: Optional[SideloadSimple],
                                   cds_markers: Optional[List[str]] = None,
                                   cds_marker_padding: int = 20000) -> SideloadedResults:
    """ Loads generic subregion/protocluster annotations from JSON files.

        Arguments:
            annotation_file: the paths to the JSON files containing annotations
            record_id: a record id to restrict annotation generation to
            manual: an optional manually specified area
            cds_markers: a list of CDS locus tags to manually generate subregions around
            cds_marker_padding: the size to include, in nucleotides, on each side of any CDS marker

        Returns:
            a GenericAnnotations instance containing all information from the
            annotation files
    """
    subregions: List[SubRegionAnnotation] = []
    protoclusters: List[ProtoclusterAnnotation] = []
    circular_origin = len(record) if record.is_circular() else None
    for annotations_file in annotation_files:
        raw = load_validated_json(annotations_file, _SCHEMA_FILE)
        tool = Tool.from_json(raw["tool"])
        for json_record in raw["records"]:
            name = json_record["name"]
            if not record.has_name(name):
                continue
            for area in json_record.get("subregions", []):
                try:
                    subregions.append(SubRegionAnnotation.from_schema_json(area, tool, circular_origin=circular_origin))
                except ValueError as err:
                    raise AntismashInputError(f"sideloaded subregion invalid: {err}")
            for area in json_record.get("protoclusters", []):
                try:
                    loaded = ProtoclusterAnnotation.from_schema_json(area, tool, circular_origin=circular_origin)
                    protoclusters.append(loaded)
                except ValueError as err:
                    raise AntismashInputError(f"sideloaded protocluster invalid: {err}")

    tool = Tool("manual", "N/A", "command line argument", {})
    if manual and record.has_name(manual.accession):
        subregion = SubRegionAnnotation(manual.start, min(manual.end, len(record.seq)), "", tool, {},
                                        circular_origin=circular_origin)
        subregions.append(subregion)

    if cds_markers:
        missing = set()
        for name in cds_markers:
            try:
                cds = record.get_cds_by_name(name)
            except KeyError:
                missing.add(name)
                continue
            start = cds.start - cds_marker_padding
            end = cds.end + cds_marker_padding
            if record.is_circular():
                start = (start + len(record)) % len(record)
                end %= len(record)
                subregion = SubRegionAnnotation(start, end, name, tool, {}, circular_origin=len(record))
            else:
                start = max(0, start)
                end = min(end, len(record))
                subregion = SubRegionAnnotation(start, end, name, tool, {})
            subregions.append(subregion)
        if missing:
            logging.warning("Features named for sideloading are not present in %s: %s", record.id, ", ".join(missing))

    for area in itertools.chain(protoclusters, subregions):
        if not record.get_cds_features_within_location(area.build_location()):
            raise AntismashInputError(f"sideloaded area contains no complete CDS features in {record.id}: {area}")

    return SideloadedResults(record.id, subregions, protoclusters)
