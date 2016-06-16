

class Constants(object):
    TASK_ID = "kinetics_tools.tasks.gather_kinetics_h5"
    DRIVER_EXE = "python -m kineticsTools.tasks.gather_kinetics_h5 --resolved-tool-contract "
    CHUNK_KEY = "$chunk.h5_id"
    OPT_CHUNK_KEY = "kinetics_tools.tasks.gather_h5_chunk_key"


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev Kinetics HDF5 Gather",
                            "General Chunk Kinetics HDF5 Gather",
                            Constants.DRIVER,
                            is_distributed=True)
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
    out = h5py.File(output_file, "w")
    first = h5py.File(chunked_files[0])
    dataLength = len(first[first.keys()[0]])
    chunkSize = min(dataLength, 8192 * 2)
    datasets = {}
    for key in first.keys():
        ds = grp.create_dataset(key, (dataLength,),
                                dtype=first[key].dtype,
                                compression="gzip",
                                chunks=(chunkSize,),
                                compression_opts=2)
        datasets[key] = ds
    for file_name in chunked_files:
        f = f5py.File(file_name)
        mask = f['base'].__array__() == ''
        for key in datasets.keys():
            datasets[key].__array__()[mask] = f[key].__array__()[mask]
    out.close()
    return f


def _run_main(chunk_input_json, output_file, chunk_key):
    chunks = load_pipeline_chunks_from_json(chunk_input_json)
    chunked_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    return gather_kinetics_h5(chunked_files, output_file)


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
