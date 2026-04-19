import argparse
from seqtools import filter_fastq


parser = argparse.ArgumentParser(
    prog="filter_fastq",
    description="Filters FASTQ file by GC%, length and mean quality. Writes output to 'filtered' subdirectory",
    epilog="This function does not return any value. The filtered reads will be written to the output FASTQ file."
)

parser.add_argument("input_fastq", type=str,
                    help="Path to the input FASTQ file to be filtered")

parser.add_argument("-o", "--output-fastq", dest="output_fastq", type=str, default="", metavar="OUT_FASTQ",
                    help="Output FASTQ name (will be created in a 'filtered' subdirectory). Default: <input>_filtered.fastq")

parser.add_argument("--overwrite",action="store_true", # pass true if mentioned
                    help="Overwrite output file if it exists")

parser.add_argument("--gc-upper", type=float, default=None, metavar="MAX_GC",
                    help="GC upper threshold (inclusive), by default: 100")
parser.add_argument("--gc-bounds", nargs=2, type=float, default=None, metavar=("MIN_GC", "MAX_GC"),
                    help="GC bounds: min max (inclusive), by default: 0 100")

parser.add_argument("--length-upper", type=int, default=None, metavar="MAX_LEN",
                    help="Length upper threshold (inclusive), by default: 2**32")
parser.add_argument("--length-bounds", nargs=2, type=int, default=None, metavar=("MIN_LEN", "MAX_LEN"),
                    help="Length bounds: min max (inclusive), by default: 0 2**32")

parser.add_argument("-q", "--quality-threshold", dest="quality_threshold", type=float, default=0, metavar="AVG_QUAL",
                    help="Minimum average read quality (inclusive), by default: 0")


if __name__ == "__main__":

    args = parser.parse_args()

    kwargs = dict(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        overwrite=args.overwrite,
        quality_threshold=args.quality_threshold,
    )

    if args.gc_upper is not None and args.gc_bounds is not None:
        print("Use either --gc-upper or --gc-bounds, not both")
        exit(1) # 1 - error, 0 - success

    if args.gc_upper is not None:
        kwargs["gc_bounds"] = args.gc_upper
    elif args.gc_bounds is not None:
        kwargs["gc_bounds"] = tuple(args.gc_bounds)
    if args.length_upper is not None and args.length_bounds is not None:
        raise ValueError("Use either --length-upper or --length-bounds, not both")

    if args.length_upper is not None and args.length_bounds is not None:
        print("Use either --length_upper or --gc-length_bounds, not both")
        exit(1) # 1 - error, 0 - success
    if args.length_upper is not None:
        kwargs["length_bounds"] = args.length_upper
    elif args.length_bounds is not None:
        kwargs["length_bounds"] = tuple(args.length_bounds)

    filter_fastq(**kwargs)