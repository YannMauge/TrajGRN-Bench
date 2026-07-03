import sys

from flecs_train import main as unified_main


if __name__ == "__main__":
    unified_main(sys.argv[1:] + ["--path_strategy", "force_timepoints"])
