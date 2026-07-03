import argparse
import sys

from flecs_train_lib import run_training


STRATEGY_CONFIG = {
    "knn": {"batch_divisor": 2},
    "force_timepoints": {"batch_divisor": 10},
}


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="inputfile", required=True)
    parser.add_argument("-o", "--outputfolder", required=True)
    parser.add_argument("-u", "--output_mode", default="full_test")
    parser.add_argument("-k", "--ko_genes", default="none")
    parser.add_argument(
        "--path_strategy",
        choices=["knn", "force_timepoints"],
        default="knn",
        help="Trajectory construction strategy.",
    )
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    config = STRATEGY_CONFIG[args.path_strategy]
    run_training(
        inputfile=args.inputfile,
        outputfolder=args.outputfolder,
        output_mode=args.output_mode,
        ko_output_genes=args.ko_genes,
        path_strategy=args.path_strategy,
        batch_divisor=config["batch_divisor"],
    )


if __name__ == "__main__":
    main(sys.argv[1:])
