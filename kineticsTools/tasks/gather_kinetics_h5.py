
"""
pbsmrtpipe task to gather chunked hdf5 files containing raw kinetics data.
"""


from collections import defaultdict
import logging
import sys

import h5py

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)


class Constants(object):
    TASK_ID = "kinetics_tools.tasks.gather_kinetics_h5"
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m kineticsTools.tasks.gather_kinetics_h5 --resolved-tool-contract "
    CHUNK_KEY = "$chunk.h5_id"
    OPT_CHUNK_KEY = "kinetics_tools.tasks.gather_h5_chunk_key"


def get_parser():
    p = get_gather_pbparser(Constants.TASK_ID,
                            Constants.VERSION,
                            "Dev Kinetics HDF5 Gather",
                            "General Chunk Kinetics HDF5 Gather",
                            Constants.DRIVER_EXE,
                            is_distributed=True,
                            default_level=logging.INFO)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with BigWig chunk key")

    p.add_output_file_type(FileTypes.H5, "h5_out",
                           "Kinetics HDF5 file",
                           "Gathered HDF5 file", "gathered")

    # Only need to add to argparse layer for the commandline
    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         Constants.CHUNK_KEY,
                         "Chunk key",
                         "Chunk key to use (format $chunk.{chunk-key}")

    return p


def gather_kinetics_h5(chunked_files, output_file):
    """
    Gather a set of 'flat' kinetics hdf5 chunks.  This will not scale well for
    large (i.e. non-bacterial) genomes.
    """
    log.info("creating {f}...".format(f=output_file))
    out = h5py.File(output_file, "w")
    first = h5py.File(chunked_files[0])
    dataLength = len(first[first.keys()[0]])
    chunkSize = min(dataLength, 8192 * 8)
    datasets = {}
    def _add_chunk(chunk):
        # FIXME this is insanely inefficient
        log.debug("  getting mask")
        mask = chunk['base'].__array__() != ''
        for key in datasets.keys():
            log.debug("  copying '{k}'".format(k=key))
            datasets[key][mask] = chunk[key][mask]
        del mask
    for key in first.keys():
        ds = out.create_dataset(key, (dataLength,),
                                dtype=first[key].dtype,
                                compression="gzip",
                                chunks=(chunkSize,),
                                compression_opts=2)
        datasets[key] = ds
    log.info("adding first chunk in {f}".format(f=chunked_files[0]))
    _add_chunk(first)
    out.flush()
    for file_name in chunked_files[1:]:
        #out = h5py.File(output_file, "r+")
        chunk = h5py.File(file_name)
        log.info("adding chunk in {f}".format(f=file_name))
        _add_chunk(chunk)
    out.close()
    return 0


def gather_kinetics_h5_byref(chunked_files, output_file):
    """
    Gather a set of hdf5 files containing with per-reference datasets.  This
    implementation sacrifices some computational efficiency to limit the
    memory overhead.
    """
    log.info("creating {f}...".format(f=output_file))
    out = h5py.File(output_file, "w")
    first = h5py.File(chunked_files[0])
    refs_by_file = defaultdict(list)
    ds_types = {}
    ref_sizes = {}
    for file_name in chunked_files:
        chunk = h5py.File(file_name)
        for rname in chunk.keys():
            refs_by_file[rname].append(file_name)
            grp = chunk[rname]
            for ds_id in grp.keys():
                if not ds_id in ds_types:
                    ds_types[ds_id] = grp[ds_id].dtype
                else:
                    assert ds_types[ds_id] == grp[ds_id].dtype
                if not rname in ref_sizes:
                    ref_sizes[rname] = len(grp[ds_id])
                else:
                    assert ref_sizes[rname] == len(grp[ds_id])
        chunk.close()
    for ref_name, file_names in refs_by_file.iteritems():
        log.info("creating group {r} (dataset length {s})".format(
                 r=ref_name,
                 s=ref_sizes[ref_name]))
        grp = out.create_group(ref_name)
        dataLength = ref_sizes[ref_name]
        chunkSize = min(dataLength, 8192)
        for ds_id, ds_type in ds_types.iteritems():
            log.debug("  dataset {i} ({t})".format(i=ds_id, t=ds_type))
            ds = grp.create_dataset(ds_id, (dataLength,), dtype=ds_type,
                                    compression="gzip", chunks=(chunkSize,))
            for file_name in file_names:
                log.debug("    reading from {f}".format(f=file_name))
                chunk = h5py.File(file_name)
                grp_chk = chunk[ref_name]
                mask = grp_chk['base'].__array__() != ''
                ds[mask] = grp_chk[ds_id][mask]
                del mask
                chunk.close()
        out.flush()
    out.close()
    return 0


def _run_main(chunk_input_json, output_file, chunk_key):
    chunks = load_pipeline_chunks_from_json(chunk_input_json)
    chunked_files = []
    for chunk in chunks:
        if chunk_key in chunk.chunk_keys:
            chunked_files.append(chunk.chunk_d[chunk_key])
        else:
            raise KeyError("Unable to find chunk key '{i}' in {p}".format(i=chunk_key, p=chunk))
    return gather_kinetics_h5_byref(chunked_files, output_file)


def args_runner(args):
    return _run_main(args.cjson_in, args.bigwig_out, args.chunk_key)


def rtc_runner(rtc):
    return _run_main(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
